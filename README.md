# stringlensing

`stringlensing` - проект для моделирования гравитационного линзирования космическими струнами, генерации синтетических изображений галактик, подгонки параметров модели и оценки качества восстановления параметров для разных телескопов, фильтров и экспозиций.

Ядро вычислений написано на C++ и экспортируется в Python через `pybind11`. Python-скрипты поверх него выполняют генерацию данных, fitting pipeline, фотометрическую подготовку FITS-изображений и телескопические simulation sweeps.

## Что находится в проекте

```text
stringlensing/
  src/stringlensing.cpp        # C++/pybind11 ядро
  setup.py                     # сборка Python extension
  fit_quality_pipeline.py      # генерация, fit и анализ качества восстановления
  telescope_simulation.py      # sweep по телескопам, фильтрам и экспозициям
  eda_photometry.py            # GUI для подготовки реальных FITS-изображений
  stringlensing_utilities.py   # вспомогательные функции визуализации/модели
params_DOT.json                # параметры DOT
params_KECK.json               # параметры KECK
requirements.txt               # Python-зависимости
ENVIRONMENT.md                 # установка окружения
BUILD_CPP_EXTENSION.md         # сборка C++ extension
docs/                          # инструкции по основным скриптам
```

## Быстрый старт

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

cd stringlensing
python setup.py build_ext --inplace --force
python -c "import stringlensing; print(stringlensing.version)"
```

После сборки можно запускать pipeline:

```bash
python fit_quality_pipeline.py --mode all --output-dir ../result_demo --n-samples 10
python telescope_simulation.py --telescopes KECK --filters V --n-samples 10 --n-exposures 2 --output-root ..
python eda_photometry.py
```

Подробные инструкции:

- [Установка окружения](ENVIRONMENT.md)
- [Сборка C++ extension](BUILD_CPP_EXTENSION.md)
- [fit_quality_pipeline.py](docs/fit_quality_pipeline.md)
- [eda_photometry.py](docs/eda_photometry.md)
- [eda_spectral.py](docs/eda_spectral.md)
- [telescope_simulation.py](docs/telescope_simulation.md)
- [Реальные FITS, CSL-1, DOT/Keck и PDF-отчёты](docs/real_data_workflow.md)
- [Карта данных и методов (PDF)](docs/reports/00_data_and_methods.pdf)
- [Анализ квазаров (PDF)](docs/reports/01_quasars.pdf)
- [CSL-1 и наклонная струна (PDF)](docs/reports/02_csl1.pdf)
- [Прогноз DOT/Keck (PDF)](docs/reports/03_telescopes.pdf)
- [Атлас всех реальных фитов (PDF)](docs/reports/04_real_fit_atlas.pdf)

## Основные возможности C++ extension

После сборки модуль `stringlensing` экспортирует функции генерации изображений и fitting utilities, включая:

- `generateUniformGalaxy_new`, `generateUniformGalaxy_inplace`
- `generateSersicGalaxy_new`, `generateSersicGalaxy_inplace`
- `generateLensedUniformGalaxy_new`, `generateLensedSersicGalaxy_new`
- `generateLensedABSersicGalaxy_new`
- `convolveWithGaussianPSF`
- `add_poisson_noise_inplace`
- `chi2_poisson`
- `levenberg_marquardt_2d`

## Что не хранится в Git

В `.gitignore` исключены тяжелые и локально воспроизводимые данные:

- `result*` - результаты прогонов, fitted JSON, plots, summaries;
- `target*` и `targets*` - FITS-данные и локальные наблюдательные выборки;
- `stringlensing/venv`, `.venv`, `venv` - виртуальные окружения;
- `stringlensing/build`, `*.so`, `*.pyd`, `*.dylib` - бинарные артефакты сборки;
- Python/cache/runtime-файлы.

Исходный `task.md`, локальные статьи и таблицы также исключены; нормализованный каталог и PDF-отчёты публикуются. Это позволяет держать репозиторий легким: код, параметры и документация попадают в Git, а большие данные остаются локально.

## Текущие исследовательские результаты

Итоговые Markdown и PDF находятся в `docs/reports/`. Инвентарь доступных FITS — в `docs/data_inventory.csv`; красные смещения с первоисточниками и расстояния по закону Хаббла — в `docs/quasar_redshifts.csv` и `docs/quasar_distances_hubble.csv`. Для воспроизведения новых cutout, фитов и графиков см. [рабочий процесс](docs/real_data_workflow.md).
