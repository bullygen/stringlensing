# `telescope_simulation.py`

`telescope_simulation.py` запускает серию прогонов `fit_quality_pipeline` для разных телескопов, фильтров и времен экспозиции. Скрипт использует параметры из `params_DOT.json` и `params_KECK.json`.

Запускать из директории `stringlensing`:

```bash
cd stringlensing
../.venv/bin/python telescope_simulation.py --help
```

## Быстрый пример

```bash
../.venv/bin/python telescope_simulation.py \
  --output-root .. \
  --telescopes KECK \
  --filters V \
  --n-samples 20 \
  --t-exp-min 10 \
  --t-exp-max 1000 \
  --n-exposures 4
```

## CLI-опции

```text
--output-root OUTPUT_ROOT
--telescopes {DOT,KECK} [{DOT,KECK} ...]
--filters FILTERS [FILTERS ...]
--n-samples N_SAMPLES
--t-exp-min T_EXP_MIN
--t-exp-max T_EXP_MAX
--n-exposures N_EXPOSURES
--seed SEED
--nx NX
--ny NY
--example-plots EXAMPLE_PLOTS
--max-retries MAX_RETRIES
--methods METHODS [METHODS ...]
--max-nfev MAX_NFEV
--max-workers MAX_WORKERS
--psf-case PSF_CASE
--psf-sigma PSF_SIGMA
```

Назначение важных опций:

- `--output-root` - корень для результатов; структура будет `<output-root>/<TELESCOPE>/<FILTER>/t_exp_<...>/`;
- `--telescopes` - список телескопов, сейчас доступны `DOT` и `KECK`;
- `--filters` - фильтры из `sky_background.filters` в JSON-конфигах;
- `--n-samples` - число synthetic samples на одну точку экспозиции;
- `--t-exp-min`, `--t-exp-max`, `--n-exposures` - логарифмическая сетка экспозиций;
- `--seed` - seed генерации;
- `--nx`, `--ny` - размер изображения;
- `--example-plots` - diagnostic plots на одну точку, capped at 10;
- `--max-workers` - число parallel workers для fit;
- `--psf-case` - имя PSF case из telescope config;
- `--psf-sigma` - ручной override sigma PSF в arcsec.

## Телескопы и фильтры

Доступные telescope configs:

```text
DOT  -> params_DOT.json
KECK -> params_KECK.json
```

Фильтры берутся из:

```json
sky_background.filters
```

Если `--filters` не указан, скрипт запускает все фильтры с числовым `sky_electrons_per_second_per_pixel`.

## PSF options

Если указан `--psf-sigma`, он имеет приоритет над конфигом.

Если указан `--psf-case`, скрипт ищет этот case в `psf.cases` каждого telescope config.

Если ничего не указано, выбирается первый доступный case из приоритетного списка:

```text
site_median_seeing
combined_seeing_0p7_plus_instrumental_center
seeing_0p7_arcsec
median_observed_V_2016_2021
```

## Выходные файлы

Для каждой точки экспозиции:

```text
<output-root>/<TELESCOPE>/<FILTER>/t_exp_<...>s/
  run_metadata.json
  synthetic_manifest.csv
  fit_manifest.csv
  synthetic/
  fit_results/
  analysis/
```

Для каждого фильтра дополнительно:

```text
<output-root>/<TELESCOPE>/<FILTER>/results/
  simulation_summary.csv
  simulation_summary.json
  success_percent_vs_t_exp.png
  recovery_success_vs_t_exp.png
  redchi_success_vs_t_exp.png
```

Результаты могут быстро становиться большими, поэтому `result*`, `target*` и `targets*` исключены из Git.

## Практические примеры

Один фильтр KECK, быстрый debug:

```bash
../.venv/bin/python telescope_simulation.py \
  --output-root .. \
  --telescopes KECK \
  --filters V \
  --n-samples 5 \
  --n-exposures 2 \
  --max-workers 2
```

Сравнение DOT и KECK:

```bash
../.venv/bin/python telescope_simulation.py \
  --output-root .. \
  --telescopes DOT KECK \
  --filters V \
  --n-samples 100 \
  --t-exp-min 10 \
  --t-exp-max 1000 \
  --n-exposures 5
```

Ручной PSF override:

```bash
../.venv/bin/python telescope_simulation.py \
  --output-root .. \
  --telescopes KECK \
  --filters V \
  --psf-sigma 0.3
```
