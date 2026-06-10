# `eda_photometry.py`

`eda_photometry.py` - графический инструмент для просмотра FITS-изображений, оценки фона/шума/PSF и экспорта реальных изображений в формат, который может читать `fit_quality_pipeline.py`.

## Запуск

```bash
cd stringlensing
../.venv/bin/python eda_photometry.py
```

Нужна графическая среда, потому что интерфейс написан на `PyQt6` и использует Matplotlib Qt backend.

## Входные данные

По умолчанию GUI смотрит в:

```text
targets_photometry/
```

Это локальная директория с FITS-данными. Она исключена из Git, потому что FITS-файлы обычно большие.

Скрипт умеет открывать FITS с несколькими HDU и многомерные массивы. Для отображения выбираются 2D-срезы с числовыми данными.

## Основной workflow

1. Запустите GUI.
2. Выберите директорию с FITS-файлами или отдельный FITS-срез.
3. Проверьте изображение и выберите область crop.
4. Укажите pixel scale, параметры crop и output folder.
5. Оцените фон, gain и PSF по изображению.
6. Экспортируйте sample в dataset format.
7. Запустите `fit_quality_pipeline.py --mode fit` по экспортированной директории.

## Выходные файлы

По умолчанию output directory:

```text
stringlensing/SL_photometry_real/
```

Для каждого экспортированного изображения создаются:

- `synthetic/images/<sample_id>.npy` - crop изображения;
- `synthetic/truth/<sample_id>.json` - placeholder/metadata JSON;
- `synthetic_manifest.csv` - manifest в формате pipeline.

Manifest содержит колонки:

```text
sample_id
image_path
truth_json_path
env_nx
env_ny
env_pixel_scale
env_noise_sigma
env_psf_sigma
```

## Запуск fit по реальным данным

После экспорта GUI показывает готовую команду. Общий вид:

```bash
cd stringlensing
../.venv/bin/python fit_quality_pipeline.py \
  --mode fit \
  --output-dir SL_photometry_real \
  --nx 100 \
  --ny 100 \
  --pixel-scale 0.2
```

Используйте `--mode fit`, потому что `--mode all` перегенерирует synthetic dataset и может заменить реальные экспортированные данные.

## Важные настройки

- `pixel-scale` должен соответствовать инструменту/данным;
- `env_noise_sigma` оценивается из фона и gain;
- `env_psf_sigma` зависит от измеренного FWHM/PSF;
- crop shape должен быть одинаковым внутри одного manifest, если не включен reset/shape override в GUI.

## Типичные проблемы

### GUI не запускается из SSH

Нужен X11/Wayland display forwarding или локальная графическая сессия.

### FITS не отображается

Проверьте, что в HDU есть числовой 2D image или многомерный массив с 2D-срезами.

### Pipeline жалуется на размер изображения

Проверьте `--nx`, `--ny` и размеры crop в `synthetic_manifest.csv`.
