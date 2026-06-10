import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import stringlensing as sl
from stringlensing_utilities import *

# ------------------ параметры для тестов ------------------

# Параметры изображения
Ny, Nx = 256, 256
pixel_scale = 0.05  # arcsec/pixel

# Параметры галактики или кружочка
Ie = 100.0
re_arcsec = 1.5
r0_arcsec = 1.5
n = 1.5
x0_pix = Nx/2 + 90
y0_pix = Ny/2 + 10
q = 0.5
pos_angle_rad = np.pi/4

'''
# Параметры струны (наклонная)
tension = 1.0e-2  # в единицах c^2/G
inclination = np.pi / 180. * 89.997  # rad
pos_angle_string_rad = np.pi/5  # rad
distane_center_string = 0.0  # arcsec
RsRg = 0.5  # отношение расстояний
'''

# Параметры струны (прямая)
tension = 1.2e-6  # в единицах c^2/G
inclination = 0.0  # rad
pos_angle_string_rad = np.pi/5  # rad
distane_center_string = 0.0  # arcsec
RsRg = 0.5  # отношение расстояний


uniform_params = {
    'Ie': Ie,
    'r0_arcsec': r0_arcsec,
    'x0_pix': x0_pix,
    'y0_pix': y0_pix
}

galaxy_params = {
    'Ie': Ie,
    're_arcsec': re_arcsec,
    'n': n,
    'x0_pix': x0_pix,
    'y0_pix': y0_pix,
    'q': q,
    'pos_angle_rad': pos_angle_rad
}

string_params = {
    'tension': tension,
    'inclination': inclination,
    'pos_angle_string_rad': pos_angle_string_rad,
    'distane_center_string': distane_center_string,
    'RsRg': RsRg
}

# ------------------ тесты прямой и обратный ------------------

def test_task_0(Ny, Nx, pixel_scale, uniform_params, string_params):
    Ie = uniform_params['Ie']
    r0_arcsec = uniform_params['r0_arcsec']
    x0_pix = uniform_params['x0_pix']
    y0_pix = uniform_params['y0_pix']

    tension = string_params['tension']
    inclination = string_params['inclination']
    pos_angle_string_rad = string_params['pos_angle_string_rad']
    distane_center_string = string_params['distane_center_string']
    RsRg = string_params['RsRg']

    # 1) In-place вариант: заранее создаём массив
    img_00 = np.empty((Ny, Nx), dtype=np.float64, order='C')
    img_10 = np.empty((Ny, Nx), dtype=np.float64, order='C')

    sl.generateUniformGalaxy_inplace(
        image = img_00,
        Ie=Ie,
        r0_arcsec=r0_arcsec,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    img_01 = sl.generateUniformGalaxy_new(
        Ny, Nx,
        Ie=Ie,
        r0_arcsec=r0_arcsec,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    img_02 = sl.convolveWithGaussianPSF(
        img_01,
        fwhm_arcsec=0.8,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    sl.generateLensedUniformGalaxy_inplace(
        image = img_10,
        Ie=Ie,
        r0_arcsec=r0_arcsec,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        tension = tension,
        inclination = inclination,
        pos_angle_string_rad = pos_angle_string_rad,
        distane_center_string = distane_center_string,
        RsRg = RsRg,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    img_11 = sl.generateLensedUniformGalaxy_new(
        Ny, Nx,
        Ie=Ie,
        r0_arcsec=r0_arcsec,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        tension = tension,
        inclination = inclination,
        pos_angle_string_rad = pos_angle_string_rad,
        distane_center_string = distane_center_string,
        RsRg = RsRg,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    sl.add_poisson_noise_inplace(img_11, background=0.0, seed=42)

    img_12 = sl.convolveWithGaussianPSF(
        img_11,
        fwhm_arcsec=0.8,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    fig, ax = plt.subplots(2, 3)
    ax[0, 0].imshow(img_00, origin='lower', cmap='gray')
    ax[0, 1].imshow(img_01, origin='lower', cmap='gray')
    ax[0, 2].imshow(img_02, origin='lower', cmap='gray')

    ax[1, 0].imshow(img_10, origin='lower', cmap='gray')
    ax[1, 1].imshow(img_11, origin='lower', cmap='gray')
    ax[1, 2].imshow(img_12, origin='lower', cmap='gray')

    plt.show()

    plot_simple_image(img_11, pixel_scale = pixel_scale, saving = False, string_params = string_params)

def test_task_1(Ny, Nx, pixel_scale, galaxy_params, string_params):

    Ie = galaxy_params['Ie']
    re_arcsec = galaxy_params['re_arcsec']
    n = galaxy_params['n']
    x0_pix = galaxy_params['x0_pix']
    y0_pix = galaxy_params['y0_pix']
    q = galaxy_params['q']
    pos_angle_rad = galaxy_params['pos_angle_rad']

    tension = string_params['tension']
    inclination = string_params['inclination']
    pos_angle_string_rad = string_params['pos_angle_string_rad']
    distane_center_string = string_params['distane_center_string']
    RsRg = string_params['RsRg']

    # 1) In-place вариант: заранее создаём массив
    img_00 = np.empty((Ny, Nx), dtype=np.float64, order='C')
    img_10 = np.empty((Ny, Nx), dtype=np.float64, order='C')

    sl.generateSersicGalaxy_inplace(
        image = img_00,
        Ie=Ie,
        re_arcsec=re_arcsec,
        n=n,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        q=q,
        pos_angle_rad = pos_angle_rad,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    img_01 = sl.generateSersicGalaxy_new(
        Ny, Nx,
        Ie=Ie,
        re_arcsec=re_arcsec,
        n=n,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        q=q,
        pos_angle_rad = pos_angle_rad,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    img_02 = sl.convolveWithGaussianPSF(
        img_01,
        fwhm_arcsec=0.8,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    sl.generateLensedSersicGalaxy_inplace(
        image = img_10,
        Ie=Ie,
        re_arcsec=re_arcsec,
        n=n,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        q=q,
        pos_angle_galaxy_rad = pos_angle_rad,
        tension = tension,
        inclination = inclination,
        pos_angle_string_rad = pos_angle_string_rad,
        distane_center_string = distane_center_string,
        RsRg = RsRg,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    img_11 = sl.generateLensedSersicGalaxy_new(
        Ny, Nx,
        Ie=Ie,
        re_arcsec=re_arcsec,
        n=n,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        q=q,
        pos_angle_galaxy_rad = pos_angle_rad,
        tension = tension,
        inclination = inclination,
        pos_angle_string_rad = pos_angle_string_rad,
        distane_center_string = distane_center_string,
        RsRg = RsRg,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    sl.add_poisson_noise_inplace(img_11, background=0.0, seed=42)

    img_12 = sl.convolveWithGaussianPSF(
        img_11,
        fwhm_arcsec=0.8,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    fig, ax = plt.subplots(2, 3)
    ax[0, 0].imshow(img_00, origin='lower', cmap='gray')
    ax[0, 1].imshow(img_01, origin='lower', cmap='gray')
    ax[0, 2].imshow(img_02, origin='lower', cmap='gray')

    ax[1, 0].imshow(img_10, origin='lower', cmap='gray')
    ax[1, 1].imshow(img_11, origin='lower', cmap='gray')
    ax[1, 2].imshow(img_12, origin='lower', cmap='gray')

    plt.show()

    plot_simple_image(img_11, pixel_scale = pixel_scale, saving = False, string_params = string_params)
    plot_fitted_image(img_12, pixel_scale = pixel_scale, saving = False, galaxy_params = galaxy_params, string_params = string_params)

    chi2 = sl.chi2_poisson(img_10, img_11, min_error=1.0)
    print(f"Chi^2 between noiseless and noisy lensed images: {chi2}")
    reduced_chi2 = chi2 / (Nx * Ny)
    print(f"Reduced Chi^2: {reduced_chi2}")

def test_task_2(Ny, Nx, pixel_scale, galaxy_params, string_params, max_iters):

    Ie = galaxy_params['Ie']
    re_arcsec = galaxy_params['re_arcsec']
    n = galaxy_params['n']
    x0_pix = galaxy_params['x0_pix']
    y0_pix = galaxy_params['y0_pix']
    q = galaxy_params['q']
    pos_angle_rad = galaxy_params['pos_angle_rad']

    tension = string_params['tension']
    inclination = string_params['inclination']
    pos_angle_string_rad = string_params['pos_angle_string_rad']
    distane_center_string = string_params['distane_center_string']
    RsRg = string_params['RsRg']

    # создание синтетической картинки
    img_syn = np.empty((Ny, Nx), dtype=np.float64, order='C')

    # линзирование
    sl.generateLensedSersicGalaxy_inplace(
        image = img_syn,
        Ie=Ie,
        re_arcsec=re_arcsec,
        n=n,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        q=q,
        pos_angle_galaxy_rad = pos_angle_rad,
        tension = tension,
        inclination = inclination,
        pos_angle_string_rad = pos_angle_string_rad,
        distane_center_string = distane_center_string,
        RsRg = RsRg,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    # шум
    sl.add_poisson_noise_inplace(img_syn, background=0.0, seed=42)

    # PSF
    img_syn_psf = sl.convolveWithGaussianPSF(
        img_syn,
        fwhm_arcsec=0.8,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    P1 = 8 * np.pi * tension * np.cos(inclination)
    P2 = 8 * np.pi * tension * np.sin(inclination)
    initial_params = Ie, re_arcsec, n, x0_pix - 10, y0_pix - 10, q, pos_angle_rad, P1, P2, pos_angle_string_rad, distane_center_string, RsRg

    result = sl.levenberg_marquardt_2d(
        model_func=stringlensing_Model,
        image=img_syn_psf.flatten(),
        Ny = Ny,
        Nx = Nx,
        initial_params=initial_params,
        max_iters=max_iters
    )
    
    print("\nРезультаты оптимизации:")
    print(f"Начальные параметры: {initial_params}")
    print(f"Конечные параметры: {result.params}")
    print(f"Приведенный χ²: {(result.chi2 / (Nx * Ny)):.4f}")
    print(f"Итераций: {result.iters}")
    print(f"Сошелся: {result.converged}")

    galaxy_params_fit = {}
    string_params_fit = {}

    galaxy_params_fit['Ie'] = result.params[0]
    galaxy_params_fit['re_arcsec'] = result.params[1]
    galaxy_params_fit['n'] = result.params[2]
    galaxy_params_fit['x0_pix'] = result.params[3]
    galaxy_params_fit['y0_pix'] = result.params[4]
    galaxy_params_fit['q'] = result.params[5]
    galaxy_params_fit['pos_angle_rad'] = result.params[6]
    string_params_fit['tension'] = np.sqrt(result.params[7]**2 + result.params[8]**2) / (8 * np.pi)
    string_params_fit['inclination'] = np.arctan(result.params[8] / result.params[7])
    string_params_fit['pos_angle_string_rad'] = result.params[9]
    string_params_fit['distane_center_string'] = result.params[10]
    string_params_fit['RsRg'] = result.params[11]

    plot_fitted_image(img_syn_psf, pixel_scale = pixel_scale, saving = False, galaxy_params = galaxy_params_fit, string_params = string_params_fit)

if __name__ == "__main__":
    test_task_0(Ny, Nx, pixel_scale, uniform_params, string_params)
    test_task_1(Ny, Nx, pixel_scale, galaxy_params, string_params)
    test_task_2(Ny, Nx, pixel_scale, galaxy_params, string_params, max_iters = 100)