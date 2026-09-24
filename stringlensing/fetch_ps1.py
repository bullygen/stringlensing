"""Инвентаризация каталога и загрузка малых r/i FITS cutout из PS1/STScI."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from astropy.io import fits
from astropy.wcs import WCS

from photometry_cli import CATALOG, ROOT, COORDINATE_OVERRIDES, catalog_rows, matching_files, source_key

HOST = 'https://ps1images.stsci.edu'
# Координаты в исходной таблице J1515+3137 противоречат имени и notes_MAST_targets.txt.
# Ниже середина пары из notes_MAST_targets.txt, до проверки по первоисточнику.
OVERRIDES = COORDINATE_OVERRIDES


def get(url: str) -> bytes:
    with urlopen(url, timeout=45) as response:
        return response.read()


def present_wcs_files(name: str, input_root: Path) -> list[Path]:
    result = []
    for p in matching_files(name, input_root):
        try:
            with fits.open(p, memmap=False, output_verify='silentfix') as hdul:
                hdul.verify('silentfix')
                if any(h.data is not None and h.data.ndim == 2 and WCS(h.header, naxis=2).has_celestial
                       for h in hdul if h.data is not None):
                    result.append(p)
        except Exception:
            continue
    return result


def retrieve(name: str, ra: float, dec: float, root: Path, size: int, filters: str):
    query = urlencode({'ra':ra,'dec':dec,'filters':filters,'type':'stack'})
    table = get(f'{HOST}/cgi-bin/ps1filenames.py?{query}').decode('utf-8').splitlines()
    if len(table) < 2:
        raise ValueError('PS1 не нашёл стек для этих координат')
    header = table[0].split()
    records = [dict(zip(header,line.split())) for line in table[1:] if line.strip()]
    saved = []
    for band in filters:
        candidates = [r for r in records if r.get('filter') == band and r.get('badflag') == '0']
        if not candidates:
            continue
        info = candidates[0]
        url = f'{HOST}/cgi-bin/fitscut.cgi?' + urlencode({'ra':ra,'dec':dec,'size':size,
                                                            'format':'fits','red':info['filename']})
        dest = root / f'PS1_cutouts_{source_key(name)}' / f'{source_key(name)}_{band}_{size}px.fits'
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            raw = get(url)
            tmp = dest.with_suffix('.tmp')
            tmp.write_bytes(raw)
            try:
                with fits.open(tmp, memmap=False) as hdul:
                    hdul.verify('silentfix')
                    sci = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
                    if not WCS(sci.header,naxis=2).has_celestial:
                        raise ValueError('Ответ без небесной WCS')
                    x,y = WCS(sci.header,naxis=2).world_to_pixel_values(ra,dec)
                    if not (0 <= x < sci.data.shape[1] and 0 <= y < sci.data.shape[0]):
                        raise ValueError('Цель вне FITS cutout')
            except Exception:
                tmp.unlink(missing_ok=True)
                raise
            tmp.replace(dest)
        provenance = {'object':name,'ra_deg':ra,'dec_deg':dec,'band':band,'url':url,
                      'source_filename':info['filename'],'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}
        dest.with_suffix('.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
        saved.append(str(dest))
    if not saved:
        raise ValueError('PS1 не нашёл пригодных r/i стеков')
    return saved


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog',type=Path,default=CATALOG)
    p.add_argument('--input-root',type=Path,default=ROOT/'targets_photometry')
    p.add_argument('--inventory',type=Path,default=ROOT/'docs'/'data_inventory.csv')
    p.add_argument('--download-missing',action='store_true')
    p.add_argument('--size',type=int,default=240)
    p.add_argument('--filters',default='ri')
    args=p.parse_args()
    output=[]
    for item in catalog_rows(args.catalog):
        name=item['object']; ra=item['ra_deg']; dec=item['dec_deg']
        note=''
        if name in OVERRIDES:
            ra,dec=OVERRIDES[name]
            note='Координаты уточнены по task.md (CSL-1) либо notes_MAST_targets.txt (J1515); требуется проверка первоисточника'
        files=matching_files(name,args.input_root)
        ready=present_wcs_files(name,args.input_root)
        downloaded=[]
        if args.download_missing and not ready and item['kind'].lower().startswith('quasar'):
            try:
                downloaded=retrieve(name,ra,dec,args.input_root,args.size,args.filters)
                ready=present_wcs_files(name,args.input_root)
                print(f'[download] {name}: {len(downloaded)} PS1 cutout',flush=True)
            except Exception as exc:
                note=(note+'; ' if note else '')+f'загрузка PS1: {exc}'
                print(f'[skip] {name}: {exc}',flush=True)
            time.sleep(0.2)
        output.append({'object':name,'kind':item['kind'],'catalog_ra_deg':item['ra_deg'],
                       'catalog_dec_deg':item['dec_deg'],'used_ra_deg':ra,'used_dec_deg':dec,
                       'local_fits_count':len(files),'wcs_fits_count':len(ready),
                       'status':'WCS FITS available' if ready else ('FITS without usable celestial WCS' if files else 'missing'),
                       'wcs_fits_paths':' | '.join(str(x.relative_to(ROOT)) for x in ready),
                       'notes':note})
    # Сохраняем также неполные строки исходной таблицы: без координат их нельзя скачать.
    with args.catalog.open(encoding='utf-8-sig',newline='') as source:
        reader=csv.DictReader(source)
        if reader.fieldnames and 'object' in reader.fieldnames:
            seen={row['object'] for row in output}
            for item in reader:
                if not item.get('object') or item['object'] in seen:
                    continue
                files=matching_files(item['object'],args.input_root)
                output.append({'object':item['object'],'kind':item.get('kind',''),
                               'catalog_ra_deg':item.get('ra_deg',''),'catalog_dec_deg':item.get('dec_deg',''),
                               'used_ra_deg':'','used_dec_deg':'','local_fits_count':len(files),
                               'wcs_fits_count':0,'status':'coordinates missing in source catalog',
                               'wcs_fits_paths':'','notes':'В исходной таблице отсутствуют численные RA/Dec; поиск FITS не выполнен'})
    args.inventory.parent.mkdir(parents=True,exist_ok=True)
    with args.inventory.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(output[0]),lineterminator='\n');writer.writeheader();writer.writerows(output)
    print(f'[inventory] {args.inventory}: objects={len(output)}, with WCS={sum(r["wcs_fits_count"]>0 for r in output)}')


if __name__=='__main__':
    main()
