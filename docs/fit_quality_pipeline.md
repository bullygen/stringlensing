# `fit_quality_pipeline.py`

`fit_quality_pipeline.py` генерирует синтетические изображения галактики со струной, подгоняет параметры модели и анализирует качество восстановления.

Запускать скрипт нужно из директории `stringlensing` после сборки C++ extension:

```bash
cd stringlensing
../.venv/bin/python fit_quality_pipeline.py --help
```

## Основные режимы

```bash
../.venv/bin/python fit_quality_pipeline.py --mode all --output-dir ../result_demo --n-samples 50
```

Доступные режимы:

- `generate` - создать synthetic dataset и `synthetic_manifest.csv`;
- `fit` - прочитать существующий synthetic dataset и выполнить fit;
- `analyze` - проанализировать результаты fit;
- `all` - выполнить `generate`, `fit`, `analyze`;
- `makeimg` - создать PNG-картинки input/model/residual из сохраненных JSON fit results.

Если нужен fit по данным, экспортированным из `eda_photometry.py`, используйте `--mode fit`, а не `--mode all`.

## CLI-опции

```text
--mode {generate,fit,analyze,all,makeimg}
--output-dir OUTPUT_DIR
--n-samples N_SAMPLES
--seed SEED
--max-retries MAX_RETRIES
--methods METHODS [METHODS ...]
--max-nfev MAX_NFEV
--nx NX
--ny NY
--pixel-scale PIXEL_SCALE
--noise-sigma NOISE_SIGMA
--psf-sigma PSF_SIGMA
--example-plots EXAMPLE_PLOTS
--makeimg
```

Назначение важных опций:

- `--output-dir` - директория результата; по умолчанию `SL_result`;
- `--n-samples` - число синтетических примеров;
- `--seed` - seed генератора случайных чисел;
- `--methods` - цепочка методов `lmfit`, по умолчанию `least_squares nelder leastsq`;
- `--max-retries` - число повторов fit для одного изображения;
- `--max-nfev` - лимит вызовов функции для оптимизатора;
- `--nx`, `--ny` - размер изображения;
- `--pixel-scale` - угловой масштаб пикселя;
- `--noise-sigma` - уровень шума в единицах изображения;
- `--psf-sigma` - sigma PSF в arcsec;
- `--example-plots` - сколько diagnostic plots сохранить;
- `--makeimg` - дополнительно сохранить картинки input/model/residual.

## Примеры

Минимальный быстрый прогон:

```bash
../.venv/bin/python fit_quality_pipeline.py \
  --mode all \
  --output-dir ../result_quick \
  --n-samples 5 \
  --max-nfev 3000
```

Только генерация:

```bash
../.venv/bin/python fit_quality_pipeline.py \
  --mode generate \
  --output-dir ../result_generated \
  --n-samples 100 \
  --seed 42
```

Fit уже созданного dataset:

```bash
../.venv/bin/python fit_quality_pipeline.py \
  --mode fit \
  --output-dir ../result_generated \
  --methods least_squares nelder leastsq \
  --max-retries 7
```

Анализ и генерация картинок:

```bash
../.venv/bin/python fit_quality_pipeline.py \
  --mode analyze \
  --output-dir ../result_generated \
  --makeimg
```

## Выходные файлы

В `--output-dir` создаются:

- `synthetic_manifest.csv` - список изображений, truth JSON и параметров окружения;
- `fit_manifest.csv` - результаты fit по каждому sample;
- `synthetic/images/*.npy` - синтетические изображения;
- `synthetic/truth/*.json` - истинные параметры модели;
- `fit_results/json/*.json` - подробные результаты fit;
- `analysis/` - summary таблицы и diagnostic plots;
- `analysis/examples/` - PNG-примеры, если включены `--example-plots` или `--makeimg`.

Директории результатов обычно тяжелые и исключены из Git правилом `result*`.

## Реальные изображения без truth

Режим `--mode real` запускает многозапусковую подгонку точечных пар из `synthetic_manifest.csv`, созданного [photometry_cli.py](real_data_workflow.md):

```bash
python fit_quality_pipeline.py --mode real \
  --real-input-dir ../result_real_prepared \
  --real-output-dir ../result_real_pointfit_v2 \
  --real-starts 6 --real-rounds 2
```

Для CSL-1 используется отдельное сравнение профилей Серсика. Модель двойного точечного источника не восстанавливает градиент разделения вдоль струны.
