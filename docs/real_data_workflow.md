# Реальные FITS, CSL-1 и PDF-отчёты

Команды выполняются из корня проекта после установки зависимостей и сборки C++ расширения по [ENVIRONMENT.md](../ENVIRONMENT.md) и [BUILD_CPP_EXTENSION.md](../BUILD_CPP_EXTENSION.md). Ниже используется готовое локальное окружение `stringlensing/venv`; в чистом checkout замените этот путь на свою `.venv/bin/python`.

```bash
cd '/home/bullygen/Documents/[Work] SAI Strings'
PY=stringlensing/venv/bin/python
$PY stringlensing/fetch_ps1.py --download-missing
$PY stringlensing/photometry_cli.py \
  --catalog docs/catalog_targets.csv \
  --size 80 --output-dir result_real_prepared \
  --exclude-object J0826+7002 \
  --exclude-object QSO1146+111B,C
$PY stringlensing/fit_quality_pipeline.py --mode real \
  --real-input-dir result_real_prepared \
  --real-output-dir result_real_pointfit_v2 \
  --real-starts 14 --real-rounds 1
```

`fetch_ps1.py` создаёт `docs/data_inventory.csv`, скачивает только объекты без пригодного локального FITS с небесной WCS, берёт PS1 r/i cutout размером 240 пикселей и кладёт рядом JSON с исходным URL и SHA256. Сервис [STScI PS1](https://outerspace.stsci.edu/spaces/PANSTARRS/pages/298812251/PS1%2BImage%2BCutout%2BService) используется без логина. Для простого обновления инвентаря запустите без `--download-missing`. Исходные большие FITS, скачанные cutout и результаты исключены из Git. Нормализованный каталог `docs/catalog_targets.csv` сохранён в Git, поэтому команды не зависят от локальной исходной таблицы.

`photometry_cli.py` принимает один FITS либо CSV-каталог. Для одного файла:

```bash
$PY stringlensing/photometry_cli.py \
  --fits targets_photometry/HST_11012006_CSL1/j9cq01010_drc.fits \
  --object CSL-1 --ra 185.877083333 --dec -12.649166667 \
  --size 128 --output-dir result_csl1_prepared
```

Для кадра с небесной WCS программа переводит RA/Dec в пиксели, измеряет WCS scale, вырезает квадрат, оценивает фон вне центральной области, переводит в e−, сохраняет `.npy`, маску плохих пикселей, JSON с WCS crop, шумом и PSF и `synthetic_manifest.csv`. Gaussian PSF получается из медианы FWHM не менее двух пригодных компактных звёздоподобных источников вне центрального объекта; `--psf-fwhm` задаёт явно измеренное внешним способом значение, и JSON укажет `user_override`. Если таких источников нет, значение PSF останется пустым: пользователь должен его измерить по другому полю/эмпирической библиотеке. В массовой предварительной подгонке тогда используются **условные** 0.10″ для HST и 1.2″ для PS1, помеченные `assumed`.

Для HST `BUNIT=ELECTRONS/S` множитель равен `EXPTIME`. Для ADU необходим gain из `CELL.GAIN`, `GAIN`, `EGAIN` или `CCDGAIN` либо аргумент `--gain`. Коэффициент `mean/variance` на уже обработанном небе использовать нельзя. Фон вычитается до сохранения; отрицательные шумовые отсчёты допустимы. Стековые PS1/HST drc изображения имеют коррелированный шум, поэтому стандартная диагональная модель ошибки приблизительна.

Кадры HCT в этой коллекции чаще имеют только линейные пиксельные координаты, а не небесную WCS. Для их включения нужны астрометрическая калибровка либо проверенные пиксельные координаты цели и масштаб. Ручной вариант: `--fits FILE --x 123 --y 456 --pixel-scale 0.3 --gain 1.5 --size 80`; значения здесь только пример, их нужно измерить для конкретного кадра. При массовом запуске без этих данных программа явно пропускает такие файлы; журнал находится в `preprocess_failures.json`. В CSV каталог CSL-1 содержит координаты, направляющие crop мимо пары: при массовом запуске применена координата из `task.md`. Координаты J1515+3137 взяты из `notes_MAST_targets.txt`, потому что CSV противоречит имени объекта; проверьте их перед научной публикацией.

## CSL-1

```bash
cd stringlensing
./venv/bin/python csl1_analysis.py \
  --input-dir ../result_csl1_prepared \
  --output-dir ../result_csl1_analysis \
  --starts 20 --fit-size 96 --gradient-max 2000
./venv/bin/python csl1_analysis.py \
  --input-dir ../result_csl1_prepared \
  --output-dir ../result_csl1_analysis_gradient2 \
  --starts 20 --fit-size 96 --gradient-max 2
./venv/bin/python csl1_analysis.py \
  --input-dir ../result_csl1_prepared \
  --output-dir ../result_csl1_analysis_gradient10 \
  --starts 40 --fit-size 96 --gradient-max 10 \
  --warm-start ../result_csl1_analysis_gradient2/summary.json
./venv/bin/python csl1_analysis.py \
  --input-dir ../result_csl1_prepared \
  --output-dir ../result_csl1_analysis \
  --starts 40 --fit-size 96 --gradient-max 2000 \
  --warm-start ../result_csl1_analysis_gradient10/summary.json
./venv/bin/python csl1_analysis.py \
  --input-dir ../result_csl1_prepared \
  --output-dir ../result_csl1_analysis_gradient2 \
  --starts 40 --fit-size 96 --gradient-max 2 \
  --warm-start ../result_csl1_analysis/summary.json
./venv/bin/python csl1_analysis.py \
  --input-dir ../result_csl1_prepared \
  --output-dir ../result_csl1_analysis_gradient10 \
  --starts 40 --fit-size 96 --gradient-max 10 \
  --warm-start ../result_csl1_analysis_gradient2/summary.json
cd ..
```

Скрипт сравнивает сдвинутую копию изображения, одну линзированную галактику Серсика и две независимые галактики Серсика. `summary.json` хранит параметры, χ² и BIC, `comparison.png` — данные/модели/вычеты. Сдвиговый тест повторяет принцип [Agol et al. (2006)](https://arxiv.org/abs/astro-ph/0603838), но не их точную маску и статистическую модель.

## DOT и Keck

```bash
$PY stringlensing/quasar_simulation.py \
  --samples 30 --magnitude 20 \
  --output-dir result_quasar_simulation_v3
```

Это новая симуляция **точечных** пар для DOT/Keck, V/R/I, двух PSF, 30/300/1000 с и разделений 0.5/1/2″. `simulation_summary.csv` содержит число успешных восстановлений, Wilson 95% CI и χ². Использованы Gaussian PSF и выбранный read noise 4 e−; это условный тест алгоритма. Старый `telescope_simulation.py` моделирует **галактики Серсика**; его отдельный повторный прогон находится в `result_telescope_rerun_20260924/`. Оба расчёта сведены в отчёте, но прогноз для точечных квазаров строится по `quasar_simulation.py`.

## Markdown и PDF

```bash
$PY docs/generate_reports.py
$PY docs/build_reports.py
```

Генератор создаёт пять `docs/reports/*.md`, PNG-графики в `docs/reports/figures/` и пять PDF рядом с Markdown. Пятый PDF `04_real_fit_atlas.pdf` содержит панель данных, моделей и вычетов **для каждого** из 44 реальных квазарных fit; для J0146–1133 нижняя панель включает переднюю галактику. В `docs/reports/data/` находятся компактные машинные сводки фитов и симуляций. PDF собираются из тех же Markdown с русским шрифтом DejaVu и формулами mathtext. Для воспроизводимости версии `markdown` и `weasyprint` указаны в `requirements.txt`.

Отчёты и таблицы должны интерпретироваться с учётом пропусков WCS/PSF и несовершенной шумовой модели. Таблица `docs/quasar_redshifts.csv` хранит ссылки на источники для 24 названных квазаров, а `docs/quasar_distances_hubble.csv` — условные расстояния по линейному закону Хаббла при H0=70 км/с/Мпк. При z порядка 1–3 это низкокрасносмещённое приближение, не космологическое расстояние. Две записи SDSS без координат остаются без идентификации. Точечная пара не измеряет градиент разделения вдоль струны; поэтому эти данные сами по себе не определяют общие $G\mu$ и $R_s$.
