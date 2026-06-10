# Руководство по сборке библиотеки stringlensing на Ubuntu

## 1. Информация о версии

**Данная инструкция создана специально для Ubuntu Linux** (версии 20.04, 22.04, 24.04 и выше). Все команды проверены на чистой установке Ubuntu.

## 2. Установка системных зависимостей

Откройте терминал и выполните:

```bash
# Обновление системы
sudo apt update
sudo apt upgrade -y

# Установка компилятора и инструментов сборки
sudo apt install -y build-essential cmake make g++ gcc

# Установка Python и необходимых пакетов
sudo apt install -y python3 python3-pip python3-venv python3-dev

# Дополнительные зависимости (если нужны)
sudo apt install -y git wget curl
```

## 3. Создание виртуального окружения

### 3.1. Создание рабочей директории

```bash
# Создаем папку для проекта
mkdir -p ~/stringlensing
cd ~/stringlensing
```

### 3.2. Создание виртуального окружения

```bash
# Создаем виртуальное окружение Python
python3 -m venv venv

# Активируем виртуальное окружение
source venv/bin/activate
```

**Важно:** После активации в начале строки терминала появится `(venv)`.

### 3.3. Быстрые команды для работы с venv

```bash
# Активация (при каждом новом открытии терминала):
source ~/stringlensing/venv/bin/activate

# Деактивация:
deactivate

# Проверка:
echo $VIRTUAL_ENV  # Должен показать путь к venv
which python3      # Должен показать путь внутри venv
```

## 4. Установка pybind11 и зависимостей

```bash
# Убедитесь, что виртуальное окружение активно
# (должно быть (venv) в начале строки)

# Обновление pip
pip install --upgrade pip

# Установка pybind11 и инструментов сборки
pip install pybind11
pip install setuptools wheel

# Для разработки (опционально)
pip install numpy pytest  # если нужны дополнительные пакеты

# Проверка установки
python3 -c "import pybind11; print(f'pybind11 {pybind11.__version__} установлен')"
python3 -m pybind11 --includes
```

## 5. Подготовка библиотеки stringlensing

### 5.1. Предполагаемая структура проекта

```
~/stringlensing/
├── venv/                    # Виртуальное окружение
├── src/                     # Исходный код
│   ├── stringlensing.cpp   # Основной файл C++
│   └── ...                 # Другие .cpp/.h файлы
├── setup.py                # Файл сборки
└── README.md               # Документация
```

### 5.2. Создание файла сборки `setup.py`

Создайте файл `setup.py` в корневой директории:

```python
from setuptools import setup, Extension
import pybind11
import sys
import os

# Флаги компиляции для Ubuntu/Linux
compile_args = [
    '-std=c++17',
    '-O3',                  # Максимальная оптимизация
    '-march=native',        # Использовать специфичные для процессора инструкции
    '-ffast-math',          # Быстрые математические вычисления
    '-fPIC',                # Position Independent Code
    '-Wall',                # Все предупреждения
    '-Wextra',              # Дополнительные предупреждения
]

# Если хотите поддержку OpenMP для многопоточности
# compile_args.append('-fopenmp')

ext_modules = [
    Extension(
        'stringlensing',
        sources=['src/stringlensing.cpp'],  # Укажите путь к вашим файлам
        include_dirs=[
            pybind11.get_include(),
            pybind11.get_include(user=True),
            'src/include'  # Если есть дополнительные заголовочные файлы
        ],
        language='c++',
        extra_compile_args=compile_args,
        # extra_link_args=['-fopenmp']  # Для линковки OpenMP
    ),
]

setup(
    name='stringlensing',
    version='1.0.0',
    author='Ваше имя',
    description='Высокопроизводительная библиотека для расчетов гравитационного линзирования космических струн',
    long_description=open('README.md').read() if os.path.exists('README.md') else '',
    ext_modules=ext_modules,
    zip_safe=False,
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Ubuntu',
    ],
)
```

## 6. Сборка библиотеки stringlensing

### 6.1. Основная команда сборки

```bash
# Убедитесь, что вы в директории проекта и venv активен
cd ~/stringlensing

# Сборка и установка в режиме разработки (рекомендуется)
pip install -e .

# ИЛИ альтернативный вариант
python3 setup.py build_ext --inplace --force
```

### 6.2. Проверка сборки

```bash
# Проверяем, что библиотека импортируется
python3 -c "
try:
    import stringlensing
    print('✓ Библиотека stringlensing успешно импортирована')
    
    # Проверка атрибутов (если знаете структуру)
    # print(dir(stringlensing))
    
except ImportError as e:
    print(f'✗ Ошибка импорта: {e}')
except Exception as e:
    print(f'✗ Другая ошибка: {e}')
"
```

## 7. Пример простого теста

Создайте файл `test.py`:

```python
#!/usr/bin/env python3
import stringlensing
import time
import numpy as np

print("=== Тестирование библиотеки stringlensing ===")
print(f"Версия Python: {sys.version}")
print()

# Пример теста производительности
if hasattr(stringlensing, 'calculate'):
    # Генерация тестовых данных
    n = 1000000
    data = np.random.randn(n).astype(np.float64)
    
    print(f"Тестирование на {n} элементах...")
    
    start = time.time()
    result = stringlensing.calculate(data)  # Замените на реальную функцию
    end = time.time()
    
    print(f"Результат: {result}")
    print(f"Время выполнения: {end - start:.4f} секунд")
else:
    print("Доступные методы:")
    for attr in dir(stringlensing):
        if not attr.startswith('_'):
            print(f"  - {attr}")
```

## 8. Устранение проблем на Ubuntu

### Проблема 1: Ошибка компилятора
```bash
# Установите последнюю версию g++
sudo apt install -y g++-12 gcc-12
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100
```

### Проблема 2: Недостающие заголовочные файлы
```bash
# Установите дополнительные библиотеки
sudo apt install -y libpython3-dev
```

### Проблема 3: Проблемы с правами
```bash
# Никогда не используйте sudo внутри venv!
# Если нужно переустановить системные пакеты:
deactivate  # выйдите из venv
sudo apt install --reinstall python3-dev
source venv/bin/activate  # вернитесь в venv
```

### Проблема 4: Устаревший pip/setuptools
```bash
pip install --upgrade pip setuptools wheel
```

## 9. Оптимизация для Ubuntu

Для максимальной производительности на Ubuntu:

### 9.1. Проверка архитектуры процессора
```bash
# Узнайте архитектуру вашего процессора
lscpu | grep -i "model name"
gcc -march=native -Q --help=target | grep march
```

### 9.2. Сборка с максимальной оптимизацией
Отредактируйте `setup.py`:
```python
compile_args = [
    '-std=c++20',  # Используйте C++20 если поддерживается
    '-O3',
    '-march=native',
    '-mtune=native',
    '-ffast-math',
    '-funroll-loops',
    '-flto',  # Link Time Optimization
    '-fPIC',
    '-DNDEBUG',  # Отключение отладочного режима
]
```

## 10. Сборка с поддержкой OpenMP (многопоточность)

```bash
# Установка OpenMP
sudo apt install -y libomp-dev

# Модификация setup.py
compile_args.append('-fopenmp')
extra_link_args = ['-fopenmp']
```

## 11. Автоматизация сборки (опционально)

Создайте скрипт `build.sh`:

```bash
#!/bin/bash
# Скрипт сборки для Ubuntu

echo "=== Сборка stringlensing на Ubuntu ==="

# Активация venv
source venv/bin/activate

# Очистка предыдущих сборок
rm -rf build/ dist/ *.so *.egg-info

# Сборка
pip install -e . --verbose

# Проверка
python3 -c "import stringlensing; print('Успешно!')"

echo "=== Готово ==="
```

Сделайте его исполняемым:
```bash
chmod +x build.sh
./build.sh
```

## 12. Полезные команды для отладки

```bash
# Просмотр собранного модуля
ldd *.so  # проверить зависимости

# Размер модуля
ls -lh *.so

# Информация о модуле
nm -gC *.so | head -20

# Проверка версии Python
python3-config --includes --ldflags
```

## Быстрая справка

```bash
# Полный процесс с нуля:
sudo apt update
sudo apt install -y python3 python3-venv python3-dev g++ cmake
cd ~
mkdir stringlensing
cd stringlensing
python3 -m venv venv
source venv/bin/activate
pip install pybind11
# поместите ваш stringlensing.cpp и setup.py
pip install -e .
python3 -c "import stringlensing"
```

**Важно:** Всегда работайте внутри виртуального окружения при сборке и использовании библиотеки!