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
    author='Igor Bulygin',
    description='A high-performance library for calculating gravitational lensing on cosmic strings',
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