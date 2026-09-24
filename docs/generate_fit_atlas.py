"""Рисунки наблюдения/модели/вычета для каждого фита реального квазара."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'docs'/'reports'
FIG = OUT/'figures'/'real_fits'


def build_atlas():
    from sys import path
    path.insert(0, str(ROOT/'stringlensing'))
    from real_fit import components, render

    with (ROOT/'result_real_prepared'/'synthetic_manifest.csv').open(newline='',encoding='utf-8') as file:
        prepared = {row['sample_id']: row for row in csv.DictReader(file)}
    with (ROOT/'result_real_pointfit_v2'/'fit_manifest.csv').open(newline='',encoding='utf-8') as file:
        fits = list(csv.DictReader(file))
    FITS_ORDER = sorted(fits, key=lambda item: (item['object'], item['sample_id']))
    FIG.mkdir(parents=True, exist_ok=True)
    report = ['# Атлас фитов реальных квазаров', '',
              f'**{len(FITS_ORDER)} изображения.** Для каждого показаны электроны после вычитания фона, '
              'модель равной пары, её вычет, маска пригодных пикселей, модель двух независимых '
              'точечных источников и её вычет. Все панели строятся из сохранённых лучших параметров '
              'фита; фиолетовый/жёлтый цвет в вычете соответствует отрицательным/положительным '
              'остаткам в единицах принятого шума. Значения χ² и BIC основаны на диагональной '
              'модели шума и поэтому являются диагностикой качества, а не точными вероятностями.', '',
              'PSF `assumed` означает условное значение: такой fit нельзя считать проверкой '
              'линзирования. Большой выигрыш независимой пары над равной часто указывает на '
              'неравные потоки, сложную PSF или смешение с другими источниками.', '']
    for index, item in enumerate(FITS_ORDER, 1):
        sid = item['sample_id']
        row = prepared[sid]
        data = np.load(row['image_path']).astype(float)
        valid = np.load(row['valid_mask_path']).astype(bool)
        result = json.loads((ROOT/'result_real_pointfit_v2'/'json'/f'{sid}.json').read_text())
        scale = float(row['env_pixel_scale'])
        fwhm = float(item['psf_fwhm_arcsec'])
        sigma = fwhm/(2*np.sqrt(2*np.log(2))*scale)
        yy, xx = np.indices(data.shape)
        equal = render(result['equal_pair']['params'], 'equal_pair', xx, yy, sigma)
        bottom_model = 'pair_galaxy' if 'pair_galaxy' in result else 'independent_pair'
        free = render(result[bottom_model]['params'], bottom_model, xx, yy, sigma)
        cx,cy=result['equal_pair']['fit_center_xy']; radius=result['equal_pair']['fit_radius_pixels']
        fit_valid=valid & ((xx-cx)**2+(yy-cy)**2 <= radius**2)
        noise = max(float(item['noise_sigma_electrons']), 1e-6)
        good = data[fit_valid]
        vmin, vmax = np.percentile(good,[2,99.7])
        if vmax <= vmin: vmax = vmin+1
        disp = lambda value: np.arcsinh(np.maximum(value-vmin,0)/(vmax-vmin)*10)
        figure, axes = plt.subplots(2,3,figsize=(9.2,6.0))
        for ax, arr, title in zip(axes.flat,
                                  [data, equal, (data-equal)/noise, fit_valid, free, (data-free)/noise],
                                  ['Данные, e−', 'Равная пара', 'Вычет равной пары / σ',
                                   'Апертура фита', 'Пара + галактика' if bottom_model=='pair_galaxy' else 'Независимая пара',
                                   'Вычет с галактикой / σ' if bottom_model=='pair_galaxy' else 'Вычет независимой / σ']):
            if 'Вычет' in title:
                ax.imshow(np.where(fit_valid,arr,np.nan),origin='lower',cmap='coolwarm',vmin=-8,vmax=8)
            elif title == 'Апертура фита':
                ax.imshow(arr,origin='lower',cmap='gray',vmin=0,vmax=1)
            else:
                ax.imshow(disp(arr),origin='lower',cmap='magma',vmin=0,vmax=np.arcsinh(10))
            ax.set_title(title,fontsize=10)
            ax.set_xticks([]);ax.set_yticks([])
        for ax, kind in [(axes[0,1],'equal_pair'),(axes[1,1],bottom_model)]:
            stars,_ = components(result[kind]['params'],kind)
            ax.scatter([p[1] for p in stars],[p[2] for p in stars],s=38,facecolors='none',
                       edgecolors='#35e7ce',linewidths=1.2)
        if bottom_model=='pair_galaxy':
            gal=result['pair_galaxy']['params']
            axes[1,1].scatter([gal[8]],[gal[9]],marker='+',s=75,c='#ffffff',linewidths=1.5)
        figure.suptitle(f"{index:02d}. {item['object']} · {Path(item['source_fits']).name}\n"
                        f"χ²ν равной пары={float(item['equal_redchi']):.2f}, независимой={float(item['free_redchi']):.2f}; "
                        f"PSF {fwhm:.2f}″ ({item['psf_status']})",fontsize=11)
        figure.subplots_adjust(left=.025,right=.975,bottom=.025,top=.87,wspace=.04,hspace=.14)
        target = FIG/f'{index:02d}_{sid}.png'
        figure.savefig(target,dpi=135)
        plt.close(figure)
        report.extend([f"## {index:02d}. {item['object']} — {Path(item['source_fits']).name}", '',
                       f"PSF: **{item['psf_status']}** ({fwhm:.3f}″); "
                       f"χ²ν равной пары **{float(item['equal_redchi']):.2f}**, независимой **{float(item['free_redchi']):.2f}**; "
                       f"ΔBIC(равная − независимая) **{float(item['delta_bic_equal_minus_free']):.1f}**; "
                       f"fit adequate: **{item['fit_adequate']}**. "
                       f"Исходный FITS: `{Path(item['source_fits']).name}`."+
                       (f" Модель с галактикой: χ²ν **{float(item['galaxy_redchi']):.2f}**, "
                        f"ΔBIC(пара − пара+галактика) **{float(item['galaxy_delta_bic_pair_minus_galaxy']):.1f}**."
                        if bottom_model=='pair_galaxy' else ''), '',
                       f'![Данные, две модели и вычеты: {item["object"]}](figures/real_fits/{target.name})', '',
                       '<div style="break-after: page"></div>', ''])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'04_real_fit_atlas.md').write_text('\n'.join(report),encoding='utf-8')
    print(f'[atlas] {len(FITS_ORDER)} figures, Markdown: {OUT/"04_real_fit_atlas.md"}')


if __name__ == '__main__':
    build_atlas()
