import stringlensing
import time
import numpy as np
import sys

print("=== Тестирование библиотеки stringlensing ===")
print(f"Версия Python: {sys.version}")
print(dir(stringlensing))

# Пример теста производительности
print("Доступные методы:")
for attr in dir(stringlensing):
    if not attr.startswith('_'):
        print(f"  - {attr}")