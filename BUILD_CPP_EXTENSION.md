# Сборка C++ extension `stringlensing`

Ядро библиотеки находится в `stringlensing/src/stringlensing.cpp` и экспортируется в Python как модуль `stringlensing` через `pybind11`.

## 1. Подготовка

Сначала активируйте окружение и установите зависимости:

```bash
cd /path/to/stringlensing-repo
source .venv/bin/activate
pip install -r requirements.txt
```

Проверьте, что `pybind11` доступен:

```bash
python -c "import pybind11; print(pybind11.__version__)"
python -m pybind11 --includes
```

## 2. Сборка на Ubuntu/Linux

```bash
cd stringlensing
../.venv/bin/python setup.py build_ext --inplace --force
```

После успешной сборки рядом с исходниками появится файл вида:

```text
stringlensing.cpython-312-x86_64-linux-gnu.so
```

Он не хранится в Git, потому что это локальный бинарный артефакт.

## 3. Проверка импорта

```bash
../.venv/bin/python -c "import stringlensing; print(stringlensing.version)"
```

Также можно посмотреть экспортированные функции:

```bash
../.venv/bin/python -c "import stringlensing; print([x for x in dir(stringlensing) if not x.startswith('_')])"
```

## 4. Альтернативная установка в editable mode

Если нужно установить extension как пакет из директории `stringlensing`:

```bash
cd stringlensing
../.venv/bin/pip install -e .
```

Для разработки чаще удобнее команда `build_ext --inplace --force`, потому что она явно пересобирает `.so` после изменений в C++.

## 5. macOS

На macOS может потребоваться Xcode Command Line Tools:

```bash
xcode-select --install
```

Затем:

```bash
cd stringlensing
../.venv/bin/python setup.py build_ext --inplace --force
```

Если стандартная сборка через `setup.py` не проходит, можно использовать ручной подход из `INSTALL_MACOS.md`, адаптировав путь к `src/stringlensing.cpp`.

## 6. Частые проблемы

### `fatal error: Python.h: No such file or directory`

Установите Python headers:

```bash
sudo apt install -y python3-dev
```

### `pybind11` не найден

```bash
source .venv/bin/activate
pip install pybind11
```

### Старый `.so` импортируется после изменений

Пересоберите с `--force`:

```bash
cd stringlensing
../.venv/bin/python setup.py build_ext --inplace --force
```

### Ошибка из-за `-march=native`

В `stringlensing/setup.py` используется `-march=native`. На некоторых машинах или при переносе бинарников между CPU это может быть неудобно. Для переносимой сборки уберите `-march=native` из `compile_args` и пересоберите extension.
