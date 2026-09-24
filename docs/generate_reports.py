"""Создание итоговых Markdown, таблиц и рисунков из результатов локальных прогонов."""
from __future__ import annotations

import csv
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.stats import binomtest

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'docs'/'reports'
FIG=OUT/'figures'
ARCSEC=206264.80624709636


def rows(path):
    with Path(path).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))


def sky_axis(fit,meta):
    header=fits.Header.fromstring(meta['wcs_header'],sep='\n')
    w=WCS(header)
    nx=int(meta['crop_origin_xy'][0]) if False else 0
    # Пиксельные координаты модели относятся к сохранённому cutout.
    data=np.load(ROOT/'result_real_prepared'/'synthetic'/'images'/f"{fit['sample_id']}.npy")
    x=(data.shape[1]-1)/2;y=(data.shape[0]-1)/2
    a=math.radians(float(fit['string_angle_image_deg']))
    p1=w.pixel_to_world(x-3*math.cos(a),y-3*math.sin(a))
    p2=w.pixel_to_world(x+3*math.cos(a),y+3*math.sin(a))
    return float(p1.position_angle(p2).rad%math.pi)


def make_figures(fits_rows,meta_by_id,csl,simulation):
    FIG.mkdir(parents=True,exist_ok=True)
    for source,name in [(ROOT/'result_csl1_analysis'/'comparison.png','csl1_comparison.png'),
                        (ROOT/'result_csl1_analysis'/'translation.png','csl1_translation.png'),
                        (ROOT/'result_csl1_analysis_gradient10'/'comparison.png','csl1_gradient10_comparison.png'),
                        (ROOT/'result_csl1_analysis_gradient2'/'comparison.png','csl1_gradient2_comparison.png'),
                        (ROOT/'result_quasar_simulation_v3'/'recovery_vs_separation.png','quasar_recovery.png')]:
        if source.exists():shutil.copy2(source,FIG/name)
    fig,ax=plt.subplots(1,2,figsize=(10,3.7))
    for status,label in [('measured','PSF по точечным источникам'),('assumed','PSF принятое')]:
        v=[float(r['equal_redchi']) for r in fits_rows if r['psf_status']==status]
        ax[0].hist(np.log10(np.clip(v,.01,1e5)),bins=np.linspace(-2,3,32),alpha=.6,label=f'{label}: {len(v)}')
    ax[0].axvline(math.log10(3),color='k',ls='--');ax[0].set_xlabel('log10 reduced chi-square');ax[0].legend(fontsize=8)
    ax[0].set_ylabel('Изображений')
    x=[float(r['delta_bic_single_minus_free']) for r in fits_rows]
    y=[float(r['delta_bic_equal_minus_free']) for r in fits_rows]
    ax[1].scatter(np.sign(x)*np.log10(1+np.abs(x)),y,c=[float(r['equal_redchi']) for r in fits_rows],s=20,cmap='viridis')
    ax[1].axhline(10,color='r',ls='--');ax[1].set_xlabel('signed log10(1+ΔBIC single−free)')
    ax[1].set_ylabel('ΔBIC equal−free');ax[1].set_ylim(-20,min(2000,max(y)+100))
    fig.tight_layout();fig.savefig(FIG/'quasar_quality.png',dpi=160);plt.close(fig)
    # Один fit на объект: PSF по точечным источникам предпочитается принятому; затем меньше redchi.
    by_object=defaultdict(list)
    for r in fits_rows:
        if r['fit_adequate']=='True' and not r['object'].startswith('QSO1146') and float(r['delta_bic_single_minus_free'])>10:
            by_object[r['object']].append(r)
    chosen={name:min(items,key=lambda r:(r['psf_status']!='measured',float(r['equal_redchi'])))
            for name,items in by_object.items()}
    angles={name:sky_axis(r,meta_by_id[r['sample_id']]) for name,r in chosen.items()}
    inventory=rows(ROOT/'docs'/'data_inventory.csv')
    fig=plt.figure(figsize=(10,5.4));ax=fig.add_subplot(111,projection='mollweide')
    for item in inventory:
        if not item['kind'].lower().startswith('quasar') or not item['used_ra_deg'] or not item['used_dec_deg']:continue
        ra=float(item['used_ra_deg']);dec=float(item['used_dec_deg'])
        x=-math.radians(((ra+180)%360)-180);y=math.radians(dec)
        ax.plot(x,y,'.',color='#64748b',ms=4)
        if item['object'] in angles:
            a=angles[item['object']];length=math.radians(5)
            dx=-length*math.sin(a)/max(math.cos(y),.15);dy=length*math.cos(a)
            ax.plot([x-dx/2,x+dx/2],[y-dy/2,y+dy/2],color='#bb3d2e',lw=1.6)
    ax.grid(alpha=.35);ax.set_title('Квазары и локальные направления пары +90°; линии не соединены')
    fig.tight_layout();fig.savefig(FIG/'quasar_mollweide.png',dpi=160);plt.close(fig)
    if angles:
        vals=np.array(list(angles.values()))
        R=float(abs(np.exp(2j*vals).mean()))
        random=np.random.default_rng(123).uniform(0,math.pi,(100000,len(vals)))
        Rnull=np.abs(np.exp(2j*random).mean(axis=1))
        p=float((1+np.count_nonzero(Rnull>=R))/(len(Rnull)+1))
    else:R=p=float('nan')
    candidates=[r for r in fits_rows if r['fit_adequate']=='True' and not r['object'].startswith('QSO1146') and
                float(r['delta_bic_single_minus_free'])>10 and float(r['delta_bic_equal_minus_free'])<10]
    redshifts={r['object']:float(r['z']) for r in rows(ROOT/'docs'/'quasar_redshifts.csv')}
    fig,ax=plt.subplots(figsize=(6.3,4.2))
    selected=sorted(candidates,key=lambda r:float(r['equal_separation_arcsec']))[:6]
    if selected:
        lim=.95*min(299792.458*redshifts[r['object']]/70/1000 for r in selected)
        rs=np.linspace(0,lim,200)
        for r in selected:
            theta=float(r['equal_separation_arcsec'])
            rg=299792.458*redshifts[r['object']]/70/1000
            ax.plot(rs,theta/(8*math.pi*ARCSEC*(1-rs/rg)),label=f"{r['object']} {theta:.2f}″")
    ax.set_yscale('log');ax.set_xlabel('$R_s$, Гпк по приближению Хаббла');ax.set_ylabel('$G\\mu$ при $i=0$')
    ax.set_title('Допустимые кривые, а не оценка общего μ')
    ax.grid(alpha=.2)
    if selected: ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(FIG/'mu_degeneracy.png',dpi=160);plt.close(fig)
    if simulation:
        fig,axgrid=plt.subplots(2,3,figsize=(10.8,4.7),sharex=True,sharey=True,constrained_layout=True)
        short_cases={'site_median_seeing':'медианный сиинг','good_site_seeing':'хороший сиинг',
                     'combined_seeing_0p7_plus_instrumental_center':'0.7″ + инструмент',
                     'instrumental_center':'инструментальная PSF'}
        for i,telescope in enumerate(['DOT','KECK']):
            for j,band in enumerate(['V','R','I']):
                ax=axgrid[i,j]
                for case in sorted({r['psf_case'] for r in simulation if r['telescope']==telescope}):
                    items=sorted((r for r in simulation if r['telescope']==telescope and r['filter']==band
                                  and r['psf_case']==case and float(r['exposure_s'])==300),
                                 key=lambda r:float(r['separation_arcsec']))
                    if not items:continue
                    x=[float(r['separation_arcsec']) for r in items]
                    y=[float(r['recovery_fraction']) for r in items]
                    lo=[y[k]-float(r['recovery_ci95_low']) for k,r in enumerate(items)]
                    hi=[float(r['recovery_ci95_high'])-y[k] for k,r in enumerate(items)]
                    ax.errorbar(x,y,yerr=[lo,hi],marker='o',capsize=2,label=short_cases.get(case,case))
                ax.set_title(f'{telescope} · {band}',fontsize=9,pad=3)
                ax.set_ylim(-.05,1.05);ax.set_xticks([.5,1,2]);ax.grid(alpha=.25)
                if i==1:ax.set_xlabel('Разделение, ″',fontsize=9)
                if j==0:ax.set_ylabel('Доля успешных фитов',fontsize=9)
                if j==2:ax.legend(fontsize=6,loc='lower right',framealpha=.8)
        fig.savefig(FIG/'quasar_recovery.png',dpi=160);plt.close(fig)
        import sys
        sys.path.insert(0,str(ROOT/'stringlensing'))
        from quasar_simulation import simulate
        rng=np.random.default_rng(20260924)
        fig,axes=plt.subplots(1,2,figsize=(8.7,4.0))
        examples=[('DOT','V','good_site_seeing',1.0),
                  ('KECK','V','combined_seeing_0p7_plus_instrumental_center',1.0)]
        for ax,(telescope,band,psf_case,separation) in zip(axes,examples):
            config=json.loads((ROOT/f'params_{telescope}.json').read_text())
            example=simulate(config,band,psf_case,300,20,separation,64,rng,return_image=True)
            data=example['example_image']
            ax.imshow(data,cmap='gray_r',origin='lower',
                      vmin=float(np.percentile(data,5)),vmax=float(np.percentile(data,99.7)),
                      interpolation='nearest')
            ax.set_title(f'{telescope} · {band} · 300 с · {separation:g}″',fontsize=9)
            ax.set_xticks([]);ax.set_yticks([])
        fig.tight_layout();fig.savefig(FIG/'quasar_sim_examples.png',dpi=150);plt.close(fig)
        fig,axs=plt.subplots(2,3,figsize=(11,4.7),sharex=True,sharey=True,constrained_layout=True)
        for i,telescope in enumerate(['DOT','KECK']):
            for j,band in enumerate(['V','R','I']):
                cases=sorted(set(r['psf_case'] for r in simulation if r['telescope']==telescope))
                case=next((c for c in cases if 'median' in c or 'combined' in c),cases[0])
                vals=[r for r in simulation if r['telescope']==telescope and r['filter']==band and r['psf_case']==case]
                tvals=sorted(set(float(r['exposure_s']) for r in vals));svals=sorted(set(float(r['separation_arcsec']) for r in vals))
                mat=np.full((len(svals),len(tvals)),np.nan)
                for r in vals:mat[svals.index(float(r['separation_arcsec'])),tvals.index(float(r['exposure_s']))]=float(r['recovery_fraction'])
                ax=axs[i,j];ax.imshow(mat,vmin=0,vmax=1,cmap='viridis',aspect='auto',origin='lower')
                ax.set_title(f'{telescope} · {band}',fontsize=9,pad=4)
                ax.set_xticks(range(len(tvals)),[f'{v:g}' for v in tvals]);ax.set_yticks(range(len(svals)),[f'{v:g}' for v in svals])
                if i==1:ax.set_xlabel('Экспозиция, с')
                if j==0:ax.set_ylabel('Разделение, ″')
        cb=fig.colorbar(axs[0,0].images[0],ax=axs.ravel().tolist(),shrink=.73,pad=.015)
        cb.set_label('Доля восстановленных пар',fontsize=9)
        fig.savefig(FIG/'quasar_sim_heatmap.png',dpi=160);plt.close(fig)
    return chosen,angles,R,p,candidates


def write_reports(inventory,prep,fitrows,csl,csl_tight,csl_physical,simulation,chosen,angles,R,p,candidates):
    OUT.mkdir(parents=True,exist_ok=True)
    n_wcs=sum(int(r['wcs_fits_count'])>0 for r in inventory)
    n_psf=sum(bool(r['env_psf_sigma']) for r in prep)
    data_table='\n'.join(f"| {r['object']} | {r['local_fits_count']} | {r['wcs_fits_count']} | {r['wcs_fits_paths'].split(' | ')[0] if r['wcs_fits_paths'] else 'нет'} |" for r in inventory)
    (OUT/'00_data_and_methods.md').write_text(f'''# Данные, код и воспроизводимость

Состояние локального исследования на 24 сентября 2026 года. Исходный каталог содержит **{len(inventory)} записей**. Для {n_wcs} объектов есть координаты и хотя бы один локальный FITS с небесной WCS; две записи SDSS1128 и SDSS111932 не содержат численных координат и остаются неразысканными. Источники: исходная таблица `spreadsheets/Stringlensing technical info - Лист1.csv`, локальные HST/HCT/PS1 и новые cutout официального [PS1 Image Cutout Service](https://outerspace.stsci.edu/spaces/PANSTARRS/pages/298812251/PS1%2BImage%2BCutout%2BService). Данные FITS и исходные статьи локальны и исключены из Git; нормализованный `docs/catalog_targets.csv`, CSV-инвентарь и инструкции входят в репозиторий.

## Карта наблюдательных данных

| Объект | FITS локально | FITS с WCS | Первый пригодный файл |
|---|---:|---:|---|
{data_table}

Полная таблица: `docs/data_inventory.csv`. В `targets_photometry/` лежат исходные снимки HST/PS1/HCT; `targets_spectra/` содержит спектры, которые в фотометрический fit не входили. Для новых PS1 cutout рядом с каждым FITS лежит JSON с URL, координатами и SHA256. Имена `J0826+7001` в HCT и `J0826+7002` в каталоге помечены как потенциальная неоднозначность. Координаты CSL-1 в CSV (RA 185.875°) не попадают на пару; использованы координаты из `task.md`: RA 185.8770833°, Dec −12.6491667°. Для J1515+3137 координаты CSV противоречат имени, использована середина пары из `notes_MAST_targets.txt` (228.915975°, +31.627875°); её нужно проверить по первоисточнику.

## Предобработка и качество

CLI `stringlensing/photometry_cli.py` выбирает SCI HDU, переводит RA/Dec в пиксели через WCS, измеряет масштаб, вырезает crop, оценивает фон на окружающем кольце и сохраняет изображения в электронах, маску плохих пикселей, PSF и полную метаинформацию. Для `ELECTRONS/S` используется EXPTIME; для ADU — gain из заголовка (`CELL.GAIN` у PS1). Если ни единицы, ни gain неизвестны, CLI отказывается выдавать фиктивные электроны. Одна оценка `mean/variance` после вычитания фона не восстанавливает gain. У стеков PS1/HST drizzle шум коррелирован, поэтому конверсия шкалы сама по себе не делает его чисто пуассоновским.

Подготовлено **{len(prep)} crop**, у **{n_psf}** Gaussian PSF оценён по двум или более соседним компактным звёздоподобным источникам; у остальных PSF не измерен и при подгонке задан явно как приближение. Ещё 26 HCT файлов не имеют небесной WCS. Для них нужны астрометрическая калибровка или проверенные пиксельные координаты и масштаб. Журнал пропусков: `result_real_prepared/preprocess_failures.json`. Оценка точности WCS проверена визуально на CSL-1: crop из CSV не содержит пару, crop по `task.md` содержит обе галактики.

## Карта ПО и команд

| Код | Назначение | Команда/инструкция |
|---|---|---|
| `stringlensing/src/stringlensing.cpp` | Профили Серсика, геометрия струны, PSF, шум | `BUILD_CPP_EXTENSION.md` |
| `stringlensing/photometry_cli.py` | FITS/WCS, фон, электроны, PSF, manifest | `docs/real_data_workflow.md` |
| `stringlensing/eda_photometry.py` | Ручной GUI-контроль crop | `docs/eda_photometry.md` |
| `stringlensing/fetch_ps1.py` | Инвентарь и загрузка PS1 cutout | `docs/real_data_workflow.md` |
| `stringlensing/fit_quality_pipeline.py` | Синтетика и точка входа `--mode real` | `docs/fit_quality_pipeline.md` |
| `stringlensing/real_fit.py` | Многозапусковая подгонка точечных пар | `docs/real_data_workflow.md` |
| `docs/generate_fit_atlas.py` | Рисунки всех реальных фитов | `docs/reports/04_real_fit_atlas.pdf` |
| `stringlensing/csl1_analysis.py` | CSL-1: сдвиговая разность и Sérsic | `docs/real_data_workflow.md` |
| `stringlensing/quasar_simulation.py` | DOT/Keck Monte Carlo точечных пар | `docs/real_data_workflow.md` |
| `stringlensing/telescope_simulation.py` | Исторический sweep галактик Серсика | `docs/telescope_simulation.md` |
| `docs/generate_reports.py`, `docs/build_reports.py` | Графики, Markdown и PDF | `docs/real_data_workflow.md` |

## Карта результатов

`result_real_prepared/` содержит cutout, маски и JSON; `result_real_pointfit_v2/` — CSV и JSON многозапусковых фитов; `result_csl1_analysis/`, `result_csl1_analysis_gradient10/` и `result_csl1_analysis_gradient2/` — сравнение моделей; `result_quasar_simulation_v3/` — новые вероятности восстановления. Старые `result_DOT/`, `result_KECK/` — прежние симуляции расширенных галактик, не подменяющие новые расчёты точечных квазаров. Публикуемые выводы и рисунки находятся в `docs/reports/` как Markdown и PDF; компактные машинные сводки каждого реального фита, всех 108 точек точечной симуляции, повторного прогона галактик и CSL-1 лежат в `docs/reports/data/`.
''',encoding='utf-8')
    meaningful_restarts=0
    for fit in fitrows:
        path=ROOT/'result_real_pointfit_v2'/'json'/f"{fit['sample_id']}.json"
        trials=json.loads(path.read_text())['equal_pair']['trial_chisqr']
        meaningful_restarts+=bool(trials and trials[0]-min(trials)>10)
    measured=sum(r['psf_status']=='measured' for r in fitrows)
    adequate=sum(r['fit_adequate']=='True' for r in fitrows)
    pair=sum(float(r['delta_bic_single_minus_free'])>10 for r in fitrows)
    good_measured=sum(r['psf_status']=='measured' for r in candidates)
    k_table='\n'.join(f"| {r['object']} | {r['psf_status']} | {float(r['equal_separation_arcsec']):.3f} | {float(r['equal_redchi']):.2f} | {float(r['delta_bic_equal_minus_free']):.1f} |" for r in candidates)
    zrows=rows(ROOT/'docs'/'quasar_redshifts.csv')
    ztable='\n'.join(f"| {r['object']} | <span style='white-space:nowrap'>{float(r['z']):.4f}</span> | <span style='white-space:nowrap'>{299792.458*float(r['z'])/70/1000:.2f}</span> | [{r['source_title']}]({r['source_url']}) | {r['match_note']} |" for r in zrows)
    (OUT/'01_quasars.md').write_text(f'''# Двойные квазары: качество подгонки и пределы вывода

Обработано **{len(fitrows)} изображения квазаров** (CSL-1 анализируется отдельно). Cutout J0826+7002 и QSO1146+111B,C исключены из всех фотометрических фитов из-за неверной области вырезки; исходные FITS и библиографические записи сохранены. Gaussian PSF оценена по звёздоподобным источникам для {measured}, для остальных {len(fitrows)-measured} использованы явно помеченные приближения: HST 0.10″, PS1 1.2″. В каждом изображении сравнивались один точечный источник, пара одинаковых PSF и пара независимых PSF; каждая пара подгонялась из 14 разных стартов с локальными ограничениями центров; для каждого кадра сохранены значения целевой функции по всем стартам и лучший найденный минимум. Для J0146–1133 отдельно подогнаны два точечных источника и экспоненциальная галактика Серсика с PSF-свёрткой. Шумовая модель учитывает shot noise и дисперсию фона, но игнорирует межпиксельные корреляции и несовершенство Gaussian PSF. Поэтому абсолютные $\\chi^2_\\nu$ и BIC служат предварительной диагностикой.

У **{meaningful_restarts}/{len(fitrows)}** изображений многозапусковый поиск улучшил первый старт по $\\chi^2$ более чем на 10; это подтверждает наличие существенных локальных минимумов, но не доказывает глобальность лучшего найденного решения.\n\nПо формальным критериям **{adequate}/{len(fitrows)}** одинаковых пар имеют $\\chi^2_\\nu<3$ и не упираются в край кадра. У **{pair}** изображений независимая пара предпочтительнее одиночного PSF на $\\Delta\\mathrm{{BIC}}>10$. Только **{len(candidates)} кадров** ({len(set(r['object'] for r in candidates))} объектов; из них {good_measured} с измеренной PSF) дополнительно имеют $\\Delta\\mathrm{{BIC}}_{{equal-free}}<10$. Это кадры для *дальнейшей проверки формы*, а не подтверждения струны.

![Распределения качества](figures/quasar_quality.png)

## Формально прошедшие численные пороги

| Объект | PSF | Разделение, ″ | $\\chi^2_\\nu$ | $\\Delta\\mathrm{{BIC}}_{{equal-free}}$ |
|---|---|---:|---:|---:|
{k_table}

Дублирующие фильтры одного объекта не являются независимыми подтверждениями. По опубликованной спектроскопии J0116+4052 является обычной галактической линзой; формальное прохождение им порога не повышает правдоподобие струнной гипотезы. Сами квазары могут иметь разную яркость; равный PSF — лишь одна проверяемая часть струнной гипотезы. На рисунке ниже показаны локальные оси, перпендикулярные разделению, без соединения разных объектов линией.

![Карта Моллвейде](figures/quasar_mollweide.png)

Для {len(angles)} объектов с обнаруженной парой и приемлемым fit осевая статистика равна $R=|\\langle e^{{2i\\phi}}\\rangle|={R:.3f}$. При моделировании 100 000 равномерных наборов получено $p={p:.3f}$. Это исследовательский тест; отбор объектов и неопределённость угла требуют отдельной коррекции. По этим данным нельзя утверждать общий протяжённый сегмент струны.

## Красные смещения и расстояния по закону Хаббла\n\n| Объект | $z$ | $D_H$, Гпк | Источник $z$ | Проверка идентификации |\n|---|---:|---:|---|---|\n{ztable}\n\nДля двух записей `SDSS1128` и `SDSS111932` из исходной таблицы нет однозначных координат и полного имени; расстояния им не приписаны. Значения $D_H$ доступны также в машинной таблице `docs/quasar_distances_hubble.csv`. Красное смещение J2015+0707 относится к координатно совпавшему объекту на расстоянии 1.7″, поэтому принадлежность обеим компонентам не подтверждена спектроскопией.\n\n## Разбор проблемных фитов\n\nУ J0116+4052, J0325-2232, J0907+0003 и J1518+4658 прежние решения с удалённым компонентом были локальными минимумами на соседних объектах. Теперь центры ограничены локальной апертурой; остатки и положения компонентов для **каждого** фильтра показаны в атласе. J0127–1441 остаётся почти неразрешённой парой: фит коллапсирует к одному PSF и не даёт надёжного разделения. У J0749+2255 и J2015+0707 потоки сильно различаются, и равная пара закономерно хуже независимой. Для J0146–1133 к паре добавлен передний профиль галактики; улучшение BIC проверяется отдельно, но остатки сохраняются. Критерии адекватности выше не превращают эти объекты в кандидатов на струну без эмпирической PSF и проверки морфологии.\n\n## Что следует из параметров модели

Для наклонной прямой струны по [Bulygin et al. 2023](https://doi.org/10.1140/epjc/s10052-023-11994-x):

$$\\theta_E=8\\pi G\\mu\\left(1-\\frac{{R_s}}{{R_g}}\\right)\\cos i$$

$$\\left.\\frac{{\\partial\\theta_E}}{{\\partial\\xi}}\\right|_0=8\\pi G\\mu\\sin i$$

Одна точечная пара определяет разделение и локальную ось, но **не определяет градиент** $\\partial\\theta_E/\\partial\\xi$: вдоль струны виден только один источник. Для всех 24 названных квазаров найдены опубликованные $z$, но у двух записей SDSS без численных координат идентификация остаётся открытой. Условная шкала расстояний вычислена **линейным законом Хаббла**: $D_H=cz/H_0$, где $c=299792.458$ км/с и $H_0=70$ км/с/Мпк. При $z\\sim1$–3 это грубая экстраполяция низкокрасносмещённого приближения; $D_H$ здесь не называется собственным или светимостным расстоянием. Подстановка $D_H$ вместо $R_g$ даёт только условную шкалу. Из одной пары всё равно нельзя получить общие $G\\mu$, $R_s$ и отдельные $i$. Даже при $i=0$ каждое измеренное разделение задаёт семейство $G\\mu=\\theta_E/[8\\pi(1-R_s/R_g)]$ (с переводом угловых секунд в радианы); при $i>0$ необходимое $G\\mu$ возрастает.

![Условные кривые Gμ при расстояниях по закону Хаббла](figures/mu_degeneracy.png)

Следующий научный шаг — получить космологически согласованные расстояния и проверяемые изображения *протяжённых* источников либо нескольких источников на одной линии, измерить PSF в каждом фильтре, после чего оценивать совместную иерархическую модель $G\\mu,R_s,i_j$.
''',encoding='utf-8')
    tr=csl['translation'];st=csl['string'];ind=csl['independent_sersic'];tight=csl_tight['string']
    physical=csl_physical['string']
    a=st['params'];b=ind['params'];c=physical['params']
    shift=math.hypot(tr['translation_dx_pix'],tr['translation_dy_pix'])
    csl_params='\n'.join([
        '| Параметр | Одна галактика + наклонная струна | Две независимые галактики |',
        '|---|---:|---:|',
        f'| Нормировка $I_e$, e на пиксель | {a[0]:.3f} | {b[0]:.3f}; {b[7]:.3f} |',
        f'| Эффективный радиус $r_e$, ″ | {a[1]:.3f} | {b[1]:.3f}; {b[8]:.3f} |',
        f'| Индекс Серсика $n$ | {a[2]:.3f} | {b[2]:.3f}; {b[9]:.3f} |',
        f'| Центр $(x_0,y_0)$, пиксель | ({a[3]:.2f}; {a[4]:.2f}) | ({b[3]:.2f}; {b[4]:.2f}); ({b[10]:.2f}; {b[11]:.2f}) |',
        f'| Отношение осей $q$ | {a[5]:.3f} | {b[5]:.3f}; {b[12]:.3f} |',
        f'| Угол профиля в координатах кадра, ° | {math.degrees(a[6]):.2f} | {math.degrees(b[6]):.2f}; {math.degrees(b[13]):.2f} |',
        f'| Разделение $\\theta_E$ при $\\xi=0$, ″ | {a[7]:.4f} | — |',
        f'| Градиент $\\partial\\theta_E/\\partial\\xi$, ″/рад | {a[8]:.3f} | — |',
        f'| Позиционный угол струны в кадре, ° | {math.degrees(a[9]):.2f} | — |',
        f'| Смещение струны, пиксель | {a[10]:.3f} | — |',
        f'| Остаточный фон, e на пиксель | {a[11]:.3f} | {b[14]:.3f} |',
        f'| $\\chi^2_\\nu$; BIC | {st["redchi"]:.3f}; {st["bic"]:.0f} | {ind["redchi"]:.3f}; {ind["bic"]:.0f} |'
    ])
    gradient_table='\n'.join([
        '| Ограничение градиента | Лучший градиент, ″/рад | $\\chi^2_\\nu$ | $\\Delta\\mathrm{BIC}$ относительно двух галактик |',
        '|---|---:|---:|---:|',
        f'| ±2000 ″/рад (свободный тест) | {a[8]:.3f} | {st["redchi"]:.3f} | {csl["delta_bic_string_minus_independent"]:.0f} |',
        f'| ±10 ″/рад | {tight["params"][8]:.3f} | {tight["redchi"]:.3f} | {csl_tight["delta_bic_string_minus_independent"]:.0f} |',
        f'| ±2 ″/рад (порядок физической границы) | {c[8]:.3f} | {physical["redchi"]:.3f} | {csl_physical["delta_bic_string_minus_independent"]:.0f} |'
    ])
    (OUT/'02_csl1.md').write_text(f'''# CSL-1: проверка модели наклонной струны

Использован HST ACS F625W `j9cq01010_drc.fits`, экспозиция 5062 с, исходная шкала `ELECTRONS/S`. По координатам из `task.md` RA 185.8770833°, Dec −12.6491667° получен WCS crop с масштабом {csl['pixel_scale_arcsec']:.4f}″/пикс, пересчитанный в электроны, и PSF FWHM {csl['psf_fwhm_arcsec']:.3f}″ по соседним звёздам. Каталожная RA 185.875° направляла crop мимо пары; этот факт проверен непосредственно изображением. В fit взято поле 96×96 пикселей вокруг обеих галактик.

## Повторение сдвигового теста

[Agol, Hogan & Plotkin (2006)](https://arxiv.org/abs/astro-ph/0603838) проверяли **прямую струну с постоянным разделением**: два изображения одной галактики должны совпасть после переноса без поворота и изменения чётности. В статье сообщено $\\chi^2_\\nu=5.8$ в F625W и 5.0 в F775W для полосы 37×100 пикселей. Здесь воспроизведена идея теста на квадратных апертурах 31×31 пиксель с оптимизацией переноса двух компонент; получилось **$\\chi^2_\\nu={tr['redchi']:.2f}$**. Число не обязано совпасть с 5.8: апертура, маска, drizzle, калибровка дисперсии и точная реализация отличаются. Разностное изображение содержит связные остатки.

Наши параметры сдвигового теста: $\\Delta x={tr['translation_dx_pix']:.3f}$ пикселя, $\\Delta y={tr['translation_dy_pix']:.3f}$ пикселя, длина сдвига **{shift:.3f} пикселя** или **{shift*csl['pixel_scale_arcsec']:.3f}″**. В статье указан сдвиг 37.9 пикселя под углом −7.1° в их повёрнутой системе координат; длины близки, но углы, апертуры и методы оценки различны. Их три свободных параметра задают проекцию струны и ширину полосы; в нашей фиксированной 31×31 апертуре оптимизируются два компонента сдвига, поэтому это воспроизведение физического теста, а не буквальной статистической процедуры.

![Сдвиговая разность](figures/csl1_translation.png)

Направления главных осей по вторым моментам в апертуре получены {csl['principal_axis_image_deg'][0]:.1f}° и {csl['principal_axis_image_deg'][1]:.1f}° в пиксельных координатах, расхождение осей **{csl['axis_difference_deg']:.1f}°**. Это согласуется с качественным выводом статьи о разных ориентациях, но наши углы нельзя напрямую приравнивать к опубликованным PA 51.7° и −2.6°: другая оценка формы и система отсчёта.

## Наклонная струна против двух галактик

### Параметры обеих моделей

Обе модели используют один HST crop, PSF FWHM {csl['psf_fwhm_arcsec']:.3f}″ и масштаб {csl['pixel_scale_arcsec']:.3f}″/пиксель. В третьем столбце параметры двух галактик разделены точкой с запятой. Координаты и углы относятся к вырезанному массиву, а не к северу на небе. Нормировка $I_e$ — значение профиля Серсика при $r_e$ в электронной шкале изображения.

{csl_params}

Оба независимых индекса $n$ упёрлись в верхнюю границу 6: сами параметры формы и абсолютная статистическая значимость неустойчивы.

### Качество подгонки

Профиль одного источника Серсика, линзированный C++ моделью с линейным $\\theta_E(\\xi)$, сравнивался с суммой двух независимых профилей Серсика при одном измеренном Gaussian PSF. По 40 разным стартам лучшая струнная модель дала $\\chi^2_\\nu={st['redchi']:.3f}$, два независимых профиля — **{ind['redchi']:.3f}**. Разница $\\mathrm{{BIC}}_{{string}}-\\mathrm{{BIC}}_{{two}}={csl['delta_bic_string_minus_independent']:.0f}$ при данной диагональной шумовой модели сильно благоприятствует двум галактикам. Лучший свободный градиент {st['params'][8]:.1f} arcsec/rad почти упёрся в предел ±2000; при ограничении ±10 arcsec/rad лучший $\\chi^2_\\nu={tight['redchi']:.3f}$. Сильное различие формы остаётся.

{gradient_table}

В [модели 2023 года](https://doi.org/10.1140/epjc/s10052-023-11994-x) градиент в радианах на радиан равен $8\\pi G\\mu\\sin i$; авторы отмечают, что при реалистичном дефицитном угле порядка $10^{{-5}}$ радиан влияние наклона на изображение мало. Свободный лучший градиент здесь соответствует {a[8]/ARCSEC:.4g} рад/рад, почти в тысячу раз больше этого масштаба. Граница ±2 ″/рад соответствует $9.7\\times10^{{-6}}$ рад/рад. Даже при таком физически более разумном ограничении два независимых профиля описывают этот кадр лучше. Формально из градиента следует $G\\mu\\geq |\\partial\\theta_E/\\partial\\xi|/(8\\pi)$, если оба угла измерять в радианах. Свободный fit требует $G\\mu\\geq {abs(a[8])/(ARCSEC*8*math.pi):.2g}$; это следствие нефизично большого градиента, а не измерение натяжения. Для CSL-1 [Agol et al.](https://arxiv.org/abs/astro-ph/0603838) приводят $z=0.463\\pm0.008$: линейная шкала Хаббла при $H_0=70$ км/с/Мпк даёт $D_H\\approx{299792.458*.463/70/1000:.2f}$ Гпк. Даже подставив её вместо $R_g$, два наблюдаемых параметра $\\theta_E$ и его градиент не определяют однозначно три неизвестные $G\\mu$, $R_s$ и $i$.

![Свободная модель наклонной струны и две галактики](figures/csl1_comparison.png)

При ограничении градиента до ±10 ″/рад остаётся следующая структура вычетов:

![Струна при ограничении градиента ±10 и две галактики](figures/csl1_gradient10_comparison.png)

При ограничении до ±2 ″/рад картина также не приближается к шумовому вычету:

![Струна при физически мотивированном ограничении градиента ±2 и две галактики](figures/csl1_gradient2_comparison.png)

### Почему сдвиговый тест не исчерпывает струнную гипотезу

[Agol et al. (2006)](https://arxiv.org/abs/astro-ph/0603838) проверяли прямую струну с **постоянным** разделением двух копий. Их статистический вывод об исключении относится именно к этой нулевой гипотезе.

1. В сдвиговом тесте один перенос совмещает всю полоску изображения. В модели наклонной струны разделение зависит от координаты вдоль неё: $\\theta_E(\\xi)=\\theta_E(0)+(\\partial\\theta_E/\\partial\\xi)_0\\xi$. Поэтому реальная наклонная линза может дать ненулевой вычет при таком тесте. Само наличие вычета не исключает эту расширенную модель.
2. Авторы ссылались на разные позиционные углы главных осей как на аргумент против линзирования. Это сильный аргумент для прямой струны с постоянным сдвигом. Наклонная или изогнутая струна может менять локальную форму изофот, так что формулировка об исключении *любого* линзирования шире параметрического семейства, проверенного сдвигом.
3. Заявленная высокая сигма характеризует несогласие с моделью постоянного переноса при принятой шумовой модели. Она не является вероятностью для всех струнных геометрий и не заменяет прямое сравнение наклонной модели с моделью двух самостоятельных галактик. Именно такое сравнение дано выше — и оно тоже оказывается неблагоприятным для струнной трактовки.

Статья 2006 года **не фитировала профили Серсика**: их основной тест был непараметрической разностью после переноса. Параметрические модели Серсика выше — наше дополнительное сравнение, а не приписываемый авторам метод.

Исследование наклонной струны из [Bulygin et al. (2023)](https://doi.org/10.1140/epjc/s10052-023-11994-x) расширяет модель прямой струны, так что тест 2006 года сам по себе не является доказательством против *всех* конфигураций. Однако количественное сравнение с этой расширенной моделью на данном кадре также не поддержало линзу. Вывод статьи 2006 года не следует объявлять ошибочным: их тест корректен для явно сформулированной модели. Остатки двухпрофильного fit тоже структурны; они могут быть следами приливного взаимодействия, несовершенством Sérsic/PSF или соседними объектами. Для различения нужны двухцветные изображения, эмпирическая PSF, ковариация drizzle и моделирование приливных компонентов.

**Ограничение значимости:** $\\chi^2$ и BIC рассчитаны с диагональной дисперсией; соседние пиксели HST drc коррелированы. Формальная разница очень велика, но точные вероятности исключения без ковариационной модели не вычислены. Это лучший найденный fit из 40 стартов, включая тёплый старт из ограниченной модели, а не математическое доказательство глобального минимума.
''',encoding='utf-8')
    if simulation:
        npoints=len(simulation);nimages=sum(int(r['n']) for r in simulation)
        sample_rows=[]
        for t in ['DOT','KECK']:
            for band in ['V','R','I']:
                group=[r for r in simulation if r['telescope']==t and r['filter']==band and float(r['exposure_s'])==300 and float(r['separation_arcsec'])==1]
                for r in group:
                    sample_rows.append(f"| {t} | {band} | {r['psf_case'].replace('site_median_seeing','медианный').replace('good_site_seeing','хороший').replace('combined_seeing_0p7_plus_instrumental_center','0.7″ + инструмент').replace('instrumental_center','инструмент')} | <span style='white-space:nowrap'>{r['n_recovered']}/{r['n']}</span> | {100*float(r['recovery_fraction']):.0f}% [{100*float(r['recovery_ci95_low']):.0f}, {100*float(r['recovery_ci95_high']):.0f}] | {float(r['median_redchi']):.2f} |")
        simtable='\n'.join(sample_rows)
    else:npoints=nimages=0;simtable='| ещё нет результатов | | | | | |'
    (OUT/'03_telescopes.md').write_text(f'''# DOT и Keck: восстановление двойных квазаров

Новый Monte Carlo использует **две точечные компоненты одинаковой яркости**, Gaussian PSF, небесный фон и масштабы пикселя из `params_DOT.json`/`params_KECK.json`. Для каждого сочетания телескопа, V/R/I фильтра, двух случаев PSF, экспозиции 30/300/1000 с и разделения 0.5/1/2″ создано 30 изображений (всего **{npoints} комбинаций, {nimages} изображений**). Оба источника условно имеют звёздную величину 20 в каждом фильтре. Поток источника выведен из указанной в JSON поверхностной яркости фона:

$$F_\\star=\\frac{{B_\\mathrm{{sky}}}}{{A_\\mathrm{{pix}}}}10^{{0.4(\\mu_\\mathrm{{sky}}-m)}}t_\\mathrm{{exp}}$$

Ширина PSF учитывает дифракцию через $\\mathrm{{FWHM}}_{{eff}}^2=\\mathrm{{FWHM}}_{{seeing}}^2+(1.028\\lambda/D)^2$ (углы в радианах) и заданное качество seeing. Шум содержит пуассоновские отсчёты источника+фона и выбранный для примера read noise 4 e−. Каждое изображение подгоняется тем же `real_fit.fit_model`, что применяется к реальным квазарам. Успех восстановления: независимая пара лучше одиночного PSF на $\\Delta\\mathrm{{BIC}}>10$, равная пара имеет $\\chi^2_\\nu<3$, ошибки разделения <20% и направления <15°. Для вероятности приведён 95% интервал Уилсона; при 30 опытах он остаётся широким. Здесь восстанавливаются разделение и направление **точечной пары**. Градиент разделения вдоль струны на одной паре не измеряется, поэтому эту долю нельзя интерпретировать как вероятность восстановления $G\\mu$, $R_s$ или наклона $i$. Отдельный прогон расширенных галактик ниже проверяет подгонку полной модели.

Примеры синтетических кадров для двух выбранных комбинаций приведены ниже в спокойной чёрно-белой шкале: слева DOT/V с хорошим сиингом, справа Keck/V с суммой атмосферного сиинга 0.7″ и инструментальной PSF. Оба кадра моделируют 300 с и пару с разделением 1″. Остальные 3240 кадров представлены сводными статистиками.

![Два примера синтетических изображений в чёрно-белой шкале](figures/quasar_sim_examples.png)

Следующий график показывает, как при фиксированной экспозиции 300 с меняется доля верно восстановленных пар при увеличении углового разделения; отдельные кривые соответствуют телескопу, фильтру и принятой PSF.

![Восстановление при 300 с](figures/quasar_recovery.png)

Ниже шесть тепловых карт: верхний ряд DOT, нижний Keck; слева направо фильтры V, R, I. В каждой карте горизонтальная ось — экспозиция 30/300/1000 с, вертикальная — расстояние между компонентами 0.5/1/2″, цвет — доля восстановленных пар от 0 до 1. Для DOT на всех картах выбран медианный сиинг площадки (`site_median_seeing`), для Keck — атмосферный сиинг 0.7″ вместе с инструментальной PSF (`combined_seeing_0p7_plus_instrumental_center`). Полные результаты для обоих вариантов PSF приведены в таблице далее.

![Доля восстановленных пар по экспозиции и разделению](figures/quasar_sim_heatmap.png)

<div style="break-before: page"></div>

## Пример при 300 с и 1″

| Телескоп | Фильтр | PSF | Успех | Доля [95% ДИ] | Медиана $\\chi^2_\\nu$ |
|---|---|---|---:|---:|---:|
{simtable}

Видимость улучшается при большем разделении и более узком PSF; экспозиция повышает S/N, но не снимает геометрическое смешение изображений. Фильтр меняет и фон, и принятую светимость источника; исходных измеренных цветов/сквозной функции не хватает для прогноза конкретного квазара. Модель не содержит wings PSF, соседних источников, насыщения, потерь от атмосферы, космических лучей и вариабельности. Поэтому числа характеризуют **условную эффективность алгоритма** при заявленных допущениях, а не гарантированную вероятность наблюдения струны.

Проверочный повторный прогон `result_telescope_rerun_20260924/` содержит по 50 симуляций расширенной галактики Серсика для DOT/V и Keck/V при 10, 31.6, 100, 316 и 1000 с. При 100 с формальный успех оптимизатора составляет 49/50 для DOT и 50/50 для Keck, но восстановление параметров в пределах 20% достигается лишь 8/50 в обоих случаях. Поэтому даже успешно завершённая оптимизация не означает физически верный fit. Архивные `result_DOT/` и `result_KECK/` тоже содержат точки по 50 симуляций **расширенных галактик Серсика**. Например, архивные DOT/U и Keck/V при 100 с дают формальный optimizer success 41/50 и 45/50, но восстановление параметров в 20% там лишь 5/50 для обеих точек. Эти массивы нельзя выдавать за прогноз для точечных квазаров. В C++ после рефакторинга координат смещение осталось нормальным к струне, поэтому численная карта эквивалентна старой в пределах округления; новое точечное Monte Carlo отвечает на отдельный наблюдательный вопрос.
''',encoding='utf-8')


def main():
    inventory=rows(ROOT/'docs'/'data_inventory.csv')
    prep=rows(ROOT/'result_real_prepared'/'synthetic_manifest.csv')
    fitrows=rows(ROOT/'result_real_pointfit_v2'/'fit_manifest.csv')
    csl=json.loads((ROOT/'result_csl1_analysis'/'summary.json').read_text())
    tight=json.loads((ROOT/'result_csl1_analysis_gradient10'/'summary.json').read_text())
    physical=json.loads((ROOT/'result_csl1_analysis_gradient2'/'summary.json').read_text())
    simpath=ROOT/'result_quasar_simulation_v3'/'simulation_summary.csv'
    simulation=rows(simpath) if simpath.exists() else []
    meta_by_id={r['sample_id']:json.loads(Path(r['metadata_path']).read_text()) for r in prep}
    chosen,angles,R,p,candidates=make_figures(fitrows,meta_by_id,csl,simulation)
    write_reports(inventory,prep,fitrows,csl,tight,physical,simulation,chosen,angles,R,p,candidates)
    data_out=OUT/'data';data_out.mkdir(parents=True,exist_ok=True)
    def export_csv(items,target):
        if not items:return
        clean=[]
        for item in items:
            item=dict(item)
            for key in ('source_fits','run_dir'):
                if key in item and item[key]:
                    try:item[key]=str(Path(item[key]).relative_to(ROOT))
                    except ValueError:item[key]=Path(item[key]).name
            clean.append(item)
        with target.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=clean[0].keys(),lineterminator='\n');writer.writeheader();writer.writerows(clean)
    export_csv(fitrows,data_out/'real_fit_summary.csv')
    trials_out=[]
    for fit in fitrows:
        payload=json.loads((ROOT/'result_real_pointfit_v2'/'json'/f"{fit['sample_id']}.json").read_text())
        for model,result in payload.items():
            for attempt,chi in enumerate(result['trial_chisqr'],1):
                trials_out.append({'sample_id':fit['sample_id'],'object':fit['object'],
                                   'model':model,'attempt':attempt,'chisqr':chi,
                                   'delta_from_best':round(chi-min(result['trial_chisqr']),3)})
    export_csv(trials_out,data_out/'fit_multistart_trials.csv')
    if simulation:export_csv(simulation,data_out/'quasar_simulation_summary.csv')
    for telescope in ('DOT','KECK'):
        source=ROOT/'result_telescope_rerun_20260924'/telescope/'V'/'results'/'simulation_summary.csv'
        if source.exists():export_csv(rows(source),data_out/f'{telescope.lower()}_extended_sersic_summary.csv')
    csl_export=dict(csl)
    csl_export['source_fits']=str(Path(csl_export['source_fits']).relative_to(ROOT))
    (data_out/'csl1_summary.json').write_text(json.dumps(csl_export,ensure_ascii=False,indent=2),encoding='utf-8')
    for label,source in [('gradient10',tight),('gradient2',physical)]:
        clean=dict(source);clean['source_fits']=str(Path(clean['source_fits']).relative_to(ROOT))
        (data_out/f'csl1_{label}_summary.json').write_text(json.dumps(clean,ensure_ascii=False,indent=2),encoding='utf-8')
    from generate_fit_atlas import build_atlas
    build_atlas()
    print(f'[reports] {len(list(OUT.glob("*.md")))} Markdown, {len(list(FIG.glob("*.png")))} summary figures')


if __name__=='__main__':main()
