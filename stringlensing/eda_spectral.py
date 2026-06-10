import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def main(dir_type, dir_source):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    bar_dir = os.path.join(script_dir, '..', dir_type)
    if len(dir_source) > 1:
        filename = dir_source
        filepath = filename if os.path.isabs(filename) else os.path.join(bar_dir, filename)
    else:
        try:
            fits_files = [f for f in os.listdir(bar_dir) if f.lower().endswith('.fits')]
            if not fits_files:
                print("Ошибка: в папке bar не найдено ни одного FITS-файла.")
                return
            filename = fits_files[0]
            filepath = os.path.join(bar_dir, filename)
            print(f"Файл не указан, использую первый найденный: {filename}")
        except FileNotFoundError:
            print(f"Ошибка: папка bar не найдена по пути {bar_dir}")
            return

    try:
        with fits.open(filepath) as hdul:
            print(f"\nОткрыт файл: {filepath}")
            hdul.info()

            data = None
            for i, hdu in enumerate(hdul):
                if hdu.data is not None:
                    data = hdu.data
                    print(f"Данные взяты из HDU #{i} ({hdu.name})")
                    return data
                    break
                else:
                    return None

            if data is None:
                print("Ошибка: в файле не найдено расширений с данными.")
                return

            

    except Exception as e:
        print(f"Ошибка при обработке файла: {e}")

dir_type = 'targets_spectra'
dir_source_1 = 'HCT_14102025_A2213/tfdcrb25.ms.fits'
dir_source_2 = 'HCT_14102025_A2213/tfdcrb26.ms.fits'
dir_source_3 = 'HCT_14102025_A2213/tfdcrb27.ms.fits'
dir_source_4 = 'HCT_14102025_A2213/tfdcrb28.ms.fits'

if __name__ == "__main__":
    data_1 = main(dir_type, dir_source_1)
    data_2 = main(dir_type, dir_source_2)
    data_3 = main(dir_type, dir_source_3)
    data_4 = main(dir_type, dir_source_4)

    fig, ax = plt.subplots(2, 2)
    ax[0, 0].plot(data_1[0, 0, :])
    ax[0, 1].plot(data_2[0, 0, :])
    ax[1, 0].plot(data_3[0, 0, :])
    ax[1, 1].plot(data_4[0, 0, :])
    plt.show()


'''
plt.figure(figsize=(10, 8))
plt.imshow(np.log(data), cmap='gray', origin='lower')
plt.colorbar(label='Интенсивность')
plt.title(f'FITS: {filename}')
plt.xlabel('X (пиксели)')
plt.ylabel('Y (пиксели)')
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 8))
plt.plot(data[0, 0, :])
plt.show()"""
'''