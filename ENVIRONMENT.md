# Установка Python-окружения

Эта инструкция готовит окружение для всех Python-скриптов проекта: `fit_quality_pipeline.py`, `telescope_simulation.py`, `eda_photometry.py` и сборки C++ extension.

## 1. Системные зависимости Ubuntu

```bash
sudo apt update
sudo apt install -y build-essential g++ gcc make cmake git
sudo apt install -y python3 python3-dev python3-pip python3-venv
sudo apt install -y libgl1 libegl1 libxcb-cursor0
```

Пакеты `libgl1`, `libegl1` и `libxcb-cursor0` полезны для `PyQt6`/Matplotlib GUI. Если запускаете только batch simulation без GUI, они могут не понадобиться.

## 2. Создание virtual environment

Рекомендуемый вариант - окружение `.venv` в корне проекта:

```bash
cd /path/to/stringlensing-repo
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

Проверка:

```bash
which python
python --version
echo "$VIRTUAL_ENV"
```

`which python` должен указывать на `.venv/bin/python`.

## 3. Установка Python-зависимостей

```bash
pip install -r requirements.txt
```

Файл `requirements.txt` содержит прямые зависимости проекта с версиями, зафиксированными из текущего рабочего окружения.

## 4. Проверка скриптов

До сборки C++ extension можно проверить CLI help:

```bash
cd stringlensing
../.venv/bin/python fit_quality_pipeline.py --help
../.venv/bin/python telescope_simulation.py --help
```

Для `eda_photometry.py` нужен графический сеанс:

```bash
../.venv/bin/python eda_photometry.py
```

## 5. Следующий шаг

После установки окружения соберите C++/Python extension по инструкции [BUILD_CPP_EXTENSION.md](BUILD_CPP_EXTENSION.md).
