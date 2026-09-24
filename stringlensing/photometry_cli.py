"""Проверяемая подготовка FITS по WCS для подгонки реальных объектов."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.detection import DAOStarFinder
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / 'docs' / 'catalog_targets.csv'
# CSV CSL-1 (185.875) misses the galaxy pair. task.md specifies 12h23m30.5s.
COORDINATE_OVERRIDES = {'CSL-1': (185.87708333333333, -12.649166666666666),
                        'J1515+3137': (228.915975, 31.627875)}
MANIFEST = ['sample_id', 'image_path', 'truth_json_path', 'env_nx', 'env_ny',
            'env_pixel_scale', 'env_noise_sigma', 'env_psf_sigma', 'metadata_path', 'valid_mask_path']


def catalog_rows(path: Path = CATALOG) -> list[dict]:
    with path.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.reader(f))
    if not rows:
        return []
    if {'object', 'kind', 'ra_deg', 'dec_deg'} <= set(rows[0]):
        columns = {name: i for i, name in enumerate(rows[0])}
        data = ([row[columns['object']], row[columns['kind']],
                 row[columns['ra_deg']], row[columns['dec_deg']]]
                for row in rows[1:] if len(row) >= len(columns))
    else:
        # Исходная таблица содержит вторую строку с подзаголовками.
        data = ([row[1], row[2], row[3], row[4]]
                for row in rows[2:] if len(row) >= 5)
    out = []
    for name, kind, ra, dec in data:
        if not name or not ra or not dec:
            continue
        try:
            out.append({'object': name, 'kind': kind,
                        'ra_deg': float(ra), 'dec_deg': float(dec)})
        except ValueError:
            continue
    return out


def source_key(name: str) -> str:
    s = name.upper().replace('–', '-').replace('—', '-')
    if 'CSL-1' in s or 'CSL1' in s:
        return 'CSL1'
    match = re.search(r'J\d{4}[+-]\d{4}', s)
    if match:
        return match.group(0)
    match = re.search(r'J\d{4}', s)
    if match:
        return match.group(0)
    if 'A2213' in s:
        return 'A2213'
    if 'QSO1146' in s:
        return 'QSO1146'
    return re.sub('[^A-Z0-9]', '', s)


def matching_files(name: str, input_root: Path) -> list[Path]:
    key = source_key(name)
    files = []
    for p in input_root.rglob('*.fits'):
        if p.name.startswith('._') or any(x in p.name.lower() for x in ['_asn.', '_flt.', '_flc.']):
            continue
        folder_key = source_key(p.parent.as_posix())
        alias = (key == 'J0826+7002' and folder_key == 'J0826+7001')
        abbreviation = (len(folder_key) == 5 and key.startswith(folder_key))
        if '/Raw/' in str(p) or (folder_key != key and not alias and not abbreviation):
            continue
        if p.name.lower().endswith(('_drz.fits', '_drc.fits')):
            if p.name.lower().endswith('_drz.fits') and p.with_name(p.name.replace('_drz.', '_drc.')).exists():
                continue
        files.append(p)
    return sorted(files)


def image_hdu(hdul: fits.HDUList, ra: float | None, dec: float | None, hdu: int | None):
    choices = [hdu] if hdu is not None else range(len(hdul))
    for i in choices:
        item = hdul[i]
        if item.data is None or item.data.ndim != 2 or min(item.data.shape) < 16:
            continue
        if item.name.upper() in {'WHT', 'ERR', 'DQ', 'CTX'}:
            continue
        header = hdul[0].header.copy()
        header.update(item.header)
        try:
            wcs = WCS(item.header, naxis=2)
            if not wcs.has_celestial:
                wcs = WCS(header, naxis=2)
        except Exception:
            wcs = None
        if ra is not None and (wcs is None or not wcs.has_celestial):
            continue
        return i, np.asarray(item.data, dtype=np.float64), header, wcs
    raise ValueError('Нет 2D SCI с небесной WCS; для снимка без WCS укажите --x и --y')


def conversion(header: fits.Header, gain_override: float | None):
    unit = str(header.get('BUNIT', '')).strip().lower().replace(' ', '')
    if unit in {'electrons/s', 'electron/s', 'e-/s', 'electrons/sec'}:
        exposure = float(header.get('EXPTIME', header.get('TEXPTIME', 0)))
        if not np.isfinite(exposure) or exposure <= 0:
            raise ValueError('ELECTRONS/S требует положительное EXPTIME')
        return exposure, 'electron_rate_times_exposure'
    if unit in {'electrons', 'electron', 'e-'}:
        return 1.0, 'electrons'
    for key in ['CELL.GAIN', 'GAIN', 'EGAIN', 'CCDGAIN', 'ATODGAIN']:
        value = gain_override if gain_override is not None else header.get(key)
        try:
            gain = float(value)
        except (ValueError, TypeError):
            continue
        if np.isfinite(gain) and gain > 0:
            return gain, f'adu_times_{"override" if gain_override is not None else key}'
    raise ValueError('Неизвестны BUNIT и gain; оценка среднего/дисперсии уже вычтенного фона ненадёжна')


def estimate_psf_fwhm(data: np.ndarray, x: float, y: float, scale: float, search_size: int = 600):
    """Медиана FWHM компактных звёзд; объект в центре исключён."""
    ny, nx = data.shape
    x0, x1 = max(0, int(x-search_size/2)), min(nx, int(x+search_size/2))
    y0, y1 = max(0, int(y-search_size/2)), min(ny, int(y+search_size/2))
    patch = np.array(data[y0:y1, x0:x1], dtype=float)
    if min(patch.shape) < 35:
        raise ValueError('Недостаточно поля для оценки PSF')
    _, median, std = sigma_clipped_stats(patch[np.isfinite(patch)], sigma=3)
    patch[~np.isfinite(patch)] = median
    yy, xx = np.indices(patch.shape)
    patch[(xx+x0-x)**2 + (yy+y0-y)**2 < 32**2] = median
    sources = DAOStarFinder(fwhm=3.0, threshold=6*std)(patch-median)
    if sources is None:
        raise ValueError('Не найдены звёзды для PSF')
    widths = []
    for star in sources[:100]:
        xcol = 'xcentroid' if 'xcentroid' in sources.colnames else 'x_centroid'
        ycol = 'ycentroid' if 'ycentroid' in sources.colnames else 'y_centroid'
        sx, sy = int(round(star[xcol])), int(round(star[ycol]))
        if sx < 7 or sy < 7 or sx+7 >= patch.shape[1] or sy+7 >= patch.shape[0]:
            continue
        stamp = patch[sy-6:sy+7, sx-6:sx+7]
        if np.max(stamp) <= median+8*std:
            continue
        gy, gx = np.indices(stamp.shape)
        def resid(p):
            amp, cx, cy, sig, bg = p
            return (bg+amp*np.exp(-((gx-cx)**2+(gy-cy)**2)/(2*sig**2))-stamp).ravel()
        try:
            fit = least_squares(resid, [float(np.max(stamp)-median), 6, 6, 1.5, median],
                                bounds=([0,3,3,0.45,-np.inf], [np.inf,9,9,8,np.inf]), max_nfev=100)
            sig = fit.x[3]
            if fit.success and 0.6 < sig < 5 and np.sqrt(np.mean(fit.fun**2)) < max(3*std, 0.25*fit.x[0]):
                widths.append(sig)
        except (ValueError, FloatingPointError):
            continue
    if len(widths) < 2:
        raise ValueError(f'Только {len(widths)} пригодных звёзд; требуется >=2')
    return float(2*math.sqrt(2*math.log(2))*np.median(widths)*scale), len(widths)


def prepare(path: Path, out: Path, name: str, ra: float | None, dec: float | None,
            size: int, hdu: int | None = None, x: float | None = None, y: float | None = None,
            gain: float | None = None, psf_fwhm: float | None = None,
            pixel_scale: float | None = None):
    if size < 16:
        raise ValueError('Размер crop должен быть >=16 пикселей')
    with fits.open(path, memmap=False, output_verify='silentfix') as hdul:
        hdul.verify('silentfix')
        i, data, header, wcs = image_hdu(hdul, ra if x is None else None, dec if y is None else None, hdu)
        if x is None or y is None:
            if ra is None or dec is None:
                raise ValueError('Нужны RA/Dec или --x/--y')
            x, y = map(float, wcs.world_to_pixel_values(ra, dec))
        if not np.isfinite([x,y]).all() or x < size/2 or y < size/2 or x >= data.shape[1]-size/2 or y >= data.shape[0]-size/2:
            raise ValueError(f'Цель вне кадра или crop у края: x={x:.2f}, y={y:.2f}, shape={data.shape}')
        if wcs is None or not wcs.has_celestial:
            if pixel_scale is None or not np.isfinite(pixel_scale) or pixel_scale <= 0:
                raise ValueError('Нет небесной WCS: для --x/--y задайте проверенный --pixel-scale')
            scale = float(pixel_scale)
            wcs = None
        else:
            scales = np.asarray(proj_plane_pixel_scales(wcs.celestial), float)*3600
            if not np.isfinite(scales).all() or min(scales) <= 0 or max(scales)/min(scales) > 1.05:
                raise ValueError(f'Некорректный или анизотропный масштаб WCS: {scales}')
            scale = float(np.mean(scales))
            if pixel_scale is not None and abs(pixel_scale/scale-1) > .05:
                raise ValueError(f'Указанный масштаб {pixel_scale} расходится с WCS {scale}')
        factor, conversion_source = conversion(header, gain)
        cut = Cutout2D(data, (x,y), (size,size), wcs=wcs, mode='strict', copy=True)
        # Оценка фона на кольце вокруг объекта, без яркого центра.
        outer = Cutout2D(data, (x,y), (min(data.shape[0], size*4), min(data.shape[1], size*4)), mode='partial', fill_value=np.nan).data
        oy, ox = np.indices(outer.shape)
        ring = (ox-(outer.shape[1]-1)/2)**2+(oy-(outer.shape[0]-1)/2)**2 > (size*0.8)**2
        sky = outer[ring & np.isfinite(outer)]
        if sky.size < 100:
            raise ValueError('Недостаточно пикселей фона')
        _, bg, scatter = sigma_clipped_stats(sky, sigma=3, maxiters=5)
        valid = np.isfinite(cut.data)
        if valid.mean() < 0.99 or not np.isfinite(scatter) or scatter <= 0:
            raise ValueError('Более 1% плохих пикселей либо невалидная оценка шума')
        electrons = np.where(valid, (cut.data-bg)*factor, 0.0)
        psf_method = 'measured_stars'
        n_stars = 0
        if psf_fwhm is None:
            try:
                psf_fwhm, n_stars = estimate_psf_fwhm(data, x, y, scale)
            except ValueError as exc:
                psf_method = f'unmeasured: {exc}'
        else:
            psf_method = 'user_override'
        safe = re.sub(r'[^A-Za-z0-9_.-]+', '_', f'{name}_{path.stem}_hdu{i}').strip('_')
        image_dir = out/'synthetic'/'images'
        meta_dir = out/'metadata'
        image_dir.mkdir(parents=True, exist_ok=True)
        meta_dir.mkdir(parents=True, exist_ok=True)
        npy = image_dir/f'{safe}.npy'
        mask_file = image_dir/f'{safe}_valid.npy'
        meta = meta_dir/f'{safe}.json'
        np.save(npy, electrons.astype(np.float32))
        np.save(mask_file, valid)
        payload = {'sample_id':safe,'object':name,'source_fits':str(path.resolve()),'hdu':i,
                   'ra_deg':ra,'dec_deg':dec,'pixel_x':x,'pixel_y':y,'crop_origin_xy':list(cut.origin_original),
                   'pixel_scale_arcsec':scale,'source_bunit':header.get('BUNIT'),
                   'conversion_to_electrons':factor,'conversion_source':conversion_source,
                   'background_input_units':float(bg),'background_sigma_input_units':float(scatter),
                   'noise_sigma_electrons':float(scatter*factor),'invalid_fraction':float(1-valid.mean()),
                   'psf_fwhm_arcsec':psf_fwhm,
                   'psf_method':psf_method,'psf_stars':n_stars,
                   'wcs_header':cut.wcs.to_header().tostring(sep='\n') if cut.wcs else None,
                   'noise_caveat':'Стек/drc имеет коррелированный шум; sigma из фона, одна лишь конверсия gain не делает его пуассоновским.'}
        meta.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
        return {'sample_id':safe,'image_path':str(npy.resolve()),'truth_json_path':'',
                'env_nx':size,'env_ny':size,'env_pixel_scale':scale,
                'env_noise_sigma':float(scatter*factor), 'env_psf_sigma':psf_fwhm or '',
                'metadata_path':str(meta.resolve()),'valid_mask_path':str(mask_file.resolve())}, payload


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fits', type=Path, help='Один FITS; иначе --catalog для всей коллекции')
    p.add_argument('--catalog', type=Path, help='CSV с объектами и координатами')
    p.add_argument('--exclude-object', action='append', default=[], help='Исключить объект с непригодным cutout; можно повторять')
    p.add_argument('--input-root', type=Path, default=ROOT/'targets_photometry')
    p.add_argument('--output-dir', type=Path, default=ROOT/'result_real_prepared')
    p.add_argument('--object', default='object')
    p.add_argument('--ra', type=float)
    p.add_argument('--dec', type=float)
    p.add_argument('--x', type=float)
    p.add_argument('--y', type=float)
    p.add_argument('--size', type=int, default=64)
    p.add_argument('--hdu', type=int)
    p.add_argument('--gain', type=float)
    p.add_argument('--pixel-scale', type=float, help='arcsec/pix для ручного crop без WCS')
    p.add_argument('--psf-fwhm', type=float)
    args=p.parse_args()
    if not args.fits and not args.catalog:
        p.error('Укажите --fits или --catalog')
    jobs=[(args.object,args.ra,args.dec,args.fits)] if args.fits else [
        (row['object'],*COORDINATE_OVERRIDES.get(row['object'],(row['ra_deg'],row['dec_deg'])),file)
        for row in catalog_rows(args.catalog) if row['object'] not in args.exclude_object
        for file in matching_files(row['object'],args.input_root)]
    rows=[]; failures=[]
    for name,ra,dec,file in jobs:
        try:
            row,meta=prepare(file,args.output_dir,name,ra,dec,args.size,args.hdu,args.x,args.y,args.gain,args.psf_fwhm,args.pixel_scale)
            rows.append(row)
            print(f"[ok] {name}: {file.name} WCS=({meta['pixel_x']:.1f},{meta['pixel_y']:.1f}) "
                  f"scale={meta['pixel_scale_arcsec']:.4f} PSF={meta['psf_fwhm_arcsec']}",flush=True)
        except Exception as exc:
            failures.append({'object':name,'fits':str(file),'error':str(exc)})
            print(f'[skip] {name}: {file}: {exc}',flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir/'synthetic_manifest.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=MANIFEST);writer.writeheader();writer.writerows(rows)
    (args.output_dir/'preprocess_failures.json').write_text(json.dumps(failures,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'[summary] prepared={len(rows)} skipped={len(failures)}')
    if not rows:
        raise SystemExit(1)


if __name__=='__main__':
    main()
