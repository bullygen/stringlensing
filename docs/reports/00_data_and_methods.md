# Данные, код и воспроизводимость

Состояние локального исследования на 24 сентября 2026 года. Исходный каталог содержит **27 записей**. Для 25 объектов есть координаты и хотя бы один локальный FITS с небесной WCS; две записи SDSS1128 и SDSS111932 не содержат численных координат и остаются неразысканными. Источники: исходная таблица `spreadsheets/Stringlensing technical info - Лист1.csv`, локальные HST/HCT/PS1 и новые cutout официального [PS1 Image Cutout Service](https://outerspace.stsci.edu/spaces/PANSTARRS/pages/298812251/PS1%2BImage%2BCutout%2BService). Данные FITS и исходные статьи локальны и исключены из Git; нормализованный `docs/catalog_targets.csv`, CSV-инвентарь и инструкции входят в репозиторий.

## Карта наблюдательных данных

| Объект | FITS локально | FITS с WCS | Первый пригодный файл |
|---|---:|---:|---|
| A2213–2652 | 15 | 1 | targets_photometry/HST_23072024_A2213–2652/HST/hst_17308_b6_acs_wfc_f814w_jf56b6/hst_17308_b6_acs_wfc_f814w_jf56b6_drc.fits |
| CSL-1 | 1 | 1 | targets_photometry/HST_11012006_CSL1/j9cq01010_drc.fits |
| DESJ0245–0556 | 2 | 2 | targets_photometry/PS1_cutouts_J0245-0556/J0245-0556_i_240px.fits |
| J0116+4052 | 5 | 5 | targets_photometry/PS1_11102010_J0116+4052/rings.v3.skycell.2161.028.stk.g.unconv.fits |
| J0127–1441 | 7 | 1 | targets_photometry/HST_30102025_J0127–1441/jf5622010_drc.fits |
| J0146–1133 | 1 | 1 | targets_photometry/HST_26102023_J0146–1133/HST/jf5628010/jf5628010_drc.fits |
| J0325-2232 | 8 | 2 | targets_photometry/PS1_cutouts_J0325-2232/J0325-2232_i_240px.fits |
| J0749+2255 | 2 | 2 | targets_photometry/HST_01052020_J0749+2255/hst_skycell-p1791x05y15_wfc3_uvis_f475w_all_drc.fits |
| J0826+7002 | 13 | 2 | targets_photometry/HST_DDMM2003_J0826+7002/jf561f010_drc.fits |
| J0907+0003 | 2 | 2 | targets_photometry/PS1_cutouts_J0907+0003/J0907+0003_i_240px.fits |
| J0907+6224 | 2 | 2 | targets_photometry/PS1_cutouts_J0907+6224/J0907+6224_i_240px.fits |
| J0937+5835 | 2 | 2 | targets_photometry/PS1_cutouts_J0937+5835/J0937+5835_i_240px.fits |
| J1515+3137 | 2 | 2 | targets_photometry/PS1_cutouts_J1515+3137/J1515+3137_i_240px.fits |
| J1518+4658 | 2 | 2 | targets_photometry/PS1_cutouts_J1518+4658/J1518+4658_i_240px.fits |
| J1524+4801 | 2 | 2 | targets_photometry/PS1_cutouts_J1524+4801/J1524+4801_i_240px.fits |
| J1540+4445 | 2 | 2 | targets_photometry/PS1_cutouts_J1540+4445/J1540+4445_i_240px.fits |
| J1616+1415 | 2 | 2 | targets_photometry/PS1_cutouts_J1616+1415/J1616+1415_i_240px.fits |
| J2015+0707 | 2 | 2 | targets_photometry/PS1_cutouts_J2015+0707/J2015+0707_i_240px.fits |
| J2032–2358 | 2 | 2 | targets_photometry/PS1_cutouts_J2032-2358/J2032-2358_i_240px.fits |
| J2132+2603 | 2 | 2 | targets_photometry/PS1_cutouts_J2132+2603/J2132+2603_i_240px.fits |
| J2147-1340 | 2 | 2 | targets_photometry/PS1_cutouts_J2147-1340/J2147-1340_i_240px.fits |
| J2250+2117 | 2 | 2 | targets_photometry/PS1_cutouts_J2250+2117/J2250+2117_i_240px.fits |
| J2316+0610 | 2 | 2 | targets_photometry/PS1_cutouts_J2316+0610/J2316+0610_i_240px.fits |
| PSJ0417+3325 | 2 | 2 | targets_photometry/PS1_cutouts_J0417+3325/J0417+3325_i_240px.fits |
| QSO1146+111B,C | 2 | 2 | targets_photometry/PS1_cutouts_QSO1146/QSO1146_i_240px.fits |
| SDSS1128 | 0 | 0 | нет |
| SDSS111932 | 0 | 0 | нет |

Полная таблица: `docs/data_inventory.csv`. В `targets_photometry/` лежат исходные снимки HST/PS1/HCT; `targets_spectra/` содержит спектры, которые в фотометрический fit не входили. Для новых PS1 cutout рядом с каждым FITS лежит JSON с URL, координатами и SHA256. Имена `J0826+7001` в HCT и `J0826+7002` в каталоге помечены как потенциальная неоднозначность. Координаты CSL-1 в CSV (RA 185.875°) не попадают на пару; использованы координаты из `task.md`: RA 185.8770833°, Dec −12.6491667°. Для J1515+3137 координаты CSV противоречат имени, использована середина пары из `notes_MAST_targets.txt` (228.915975°, +31.627875°); её нужно проверить по первоисточнику.

## Предобработка и качество

CLI `stringlensing/photometry_cli.py` выбирает SCI HDU, переводит RA/Dec в пиксели через WCS, измеряет масштаб, вырезает crop, оценивает фон на окружающем кольце и сохраняет изображения в электронах, маску плохих пикселей, PSF и полную метаинформацию. Для `ELECTRONS/S` используется EXPTIME; для ADU — gain из заголовка (`CELL.GAIN` у PS1). Если ни единицы, ни gain неизвестны, CLI отказывается выдавать фиктивные электроны. Одна оценка `mean/variance` после вычитания фона не восстанавливает gain. У стеков PS1/HST drizzle шум коррелирован, поэтому конверсия шкалы сама по себе не делает его чисто пуассоновским.

Подготовлено **45 crop**, у **25** Gaussian PSF оценён по двум или более соседним компактным звёздоподобным источникам; у остальных PSF не измерен и при подгонке задан явно как приближение. Ещё 26 HCT файлов не имеют небесной WCS. Для них нужны астрометрическая калибровка или проверенные пиксельные координаты и масштаб. Журнал пропусков: `result_real_prepared/preprocess_failures.json`. Оценка точности WCS проверена визуально на CSL-1: crop из CSV не содержит пару, crop по `task.md` содержит обе галактики.

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
