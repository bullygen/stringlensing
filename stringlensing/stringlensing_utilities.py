import numpy as np
import matplotlib.pyplot as plt
import stringlensing as sl

def get_line_coords(image_shape, pixel_scale, theta, h):
    ny, nx = image_shape
    # Центр изображения в пикселях
    cx = (nx - 1) / 2.0
    cy = (ny - 1) / 2.0
    # Границы в "центральной" системе координат (x,y), где (0,0) в центре
    x_min, x_max = -cx, cx
    y_min, y_max = -cy, cy
    # Прямая: a x + b y = c
    a =  - np.sin(theta)   # компонент нормали по x
    b = np.cos(theta)   # компонент нормали по y
    c = h / pixel_scale  # расстояние от центра до прямой
    points = []
    # Пересечение с левым краем: x = x_min
    if abs(b) > 1e-12:
        y = (c - a * x_min) / b
        if y_min - 1e-9 <= y <= y_max + 1e-9:
            points.append((x_min, y))
    # Пересечение с правым краем: x = x_max
    if abs(b) > 1e-12:
        y = (c - a * x_max) / b
        if y_min - 1e-9 <= y <= y_max + 1e-9:
            points.append((x_max, y))
    # Пересечение с нижним краем: y = y_min
    if abs(a) > 1e-12:
        x = (c - b * y_min) / a
        if x_min - 1e-9 <= x <= x_max + 1e-9:
            points.append((x, y_min))
    # Пересечение с верхним краем: y = y_max
    if abs(a) > 1e-12:
        x = (c - b * y_max) / a
        if x_min - 1e-9 <= x <= x_max + 1e-9:
            points.append((x, y_max))
    '''
    # Убираем дубликаты (может пересечь в углу)
    unique_points = []
    for p in points:
        if not any(np.allclose(p, q) for q in unique_points):
            unique_points.append(p)
    '''
    
    # Берём первые две уникальные точки
    (x1, y1), (x2, y2) = points[:2]
    
    # Переводим обратно в координаты пикселей:
    # x_pix = x + cx, y_pix = y + cy
    x1_pix, y1_pix = x1 * pixel_scale, y1 * pixel_scale
    x2_pix, y2_pix = x2 * pixel_scale, y2 * pixel_scale
    return [x1_pix, x2_pix], [y1_pix, y2_pix]

def get_string_line_coords(image_shape, pixel_scale, string_params):
    string_position = {}
    theta_0 = string_params['pos_angle_string_rad']
    h_0 = string_params['distane_center_string'] * pixel_scale  # C++ задаёт смещение в пикселях
    Re = 8 * np.pi * string_params['tension'] * 206265 * (1 - string_params['RsRg']) * np.cos(string_params['inclination'])
    alpha = np.arctan(8 * np.pi * string_params['tension'] * np.sin(string_params['inclination']))
    h_plus = (h_0 + Re) * np.cos(alpha)
    h_minus = (h_0 - Re) * np.cos(alpha)
    theta_plus = theta_0 + alpha
    theta_minus = theta_0 - alpha

    string_position['x0'], string_position['y0'] = get_line_coords(
        image_shape,
        pixel_scale,
        theta_0,
        h_0
    )

    string_position['x+'], string_position['y+'] = get_line_coords(
        image_shape,
        pixel_scale,
        theta_plus,
        h_plus
    )

    string_position['x-'], string_position['y-'] = get_line_coords(
        image_shape,
        pixel_scale,
        theta_minus,
        h_minus
    )
    
    return string_position

def plot_simple_image(image, pixel_scale, saving, string_params):
    image_shape = image.shape
    
    extent = [
        -image.shape[1] / 2 * pixel_scale,
         image.shape[1] / 2 * pixel_scale,
        -image.shape[0] / 2 * pixel_scale,
         image.shape[0] / 2 * pixel_scale
    ]
    
    fig, ax = plt.subplots(1, 2, figsize=(6, 3))

    # Plot the image
    ax[0].imshow(image, origin='lower', cmap='gray', extent = extent)
    ax[0].set_xlabel('Arcseconds')
    ax[0].set_ylabel('Arcseconds')
    ax[0].set_title('Lensed Sersic Galaxy Image')

    string_position = get_string_line_coords(image_shape, pixel_scale, string_params)
    ax[1].imshow(image, origin='lower', cmap='gray', extent = extent)
    ax[1].plot(string_position['x0'], string_position['y0'], color='red')
    ax[1].plot(string_position['x+'], string_position['y+'], 'r--')
    ax[1].plot(string_position['x-'], string_position['y-'], 'r--')
    ax[1].set_xlabel('Arcseconds')
    ax[1].set_ylabel('Arcseconds')
    ax[1].set_title('Cosmic String Position')

    if saving:
        plt.imsave('lensed_sersic_galaxy.png', image)

    plt.show()
 
def plot_fitted_image(image, pixel_scale, saving, galaxy_params, string_params):
    image_shape = image.shape 
    fwhm_arcsec = 0.8

    extent = [
        -image.shape[1] / 2 * pixel_scale,
         image.shape[1] / 2 * pixel_scale,
        -image.shape[0] / 2 * pixel_scale,
         image.shape[0] / 2 * pixel_scale
    ]
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # Plot the real data image
    ax[0, 0].imshow(image, origin='lower', extent=extent, cmap='gray')
    #ax[0, 0].colorbar(label='Intensity')
    ax[0, 0].set_xlabel('Arcseconds')
    ax[0, 0].set_ylabel('Arcseconds')
    ax[0, 0].set_title('Candidate for Lensed Sersic Galaxy')

    # Plot the model (lensed)
    model_ = sl.generateLensedSersicGalaxy_new(
        image_shape[0], image_shape[1],
        Ie= galaxy_params['Ie'],
        re_arcsec= galaxy_params['re_arcsec'],
        n= galaxy_params['n'],
        x0_pix= galaxy_params['x0_pix'],
        y0_pix= galaxy_params['y0_pix'],
        q= galaxy_params['q'],
        pos_angle_galaxy_rad = galaxy_params['pos_angle_rad'],
        tension = string_params['tension'],
        inclination = string_params['inclination'],
        pos_angle_string_rad = string_params['pos_angle_string_rad'],
        distane_center_string = string_params['distane_center_string'],
        RsRg = string_params['RsRg'],
        pixel_scale_arcsec_per_pix=pixel_scale
    )
    model = sl.convolveWithGaussianPSF(
        model_,
        fwhm_arcsec=fwhm_arcsec,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    ax[0, 1].imshow(model, origin='lower', extent=extent, cmap='gray')
    #ax[0, 1].colorbar(label='Intensity')
    ax[0, 1].set_xlabel('Arcseconds')
    ax[0, 1].set_ylabel('Arcseconds')
    ax[0, 1].set_title('Modelled Lensed Sersic Galaxy')

    # Plot the real data image + the string
    string_position = get_string_line_coords(image_shape, pixel_scale, string_params)
    ax[1, 0].imshow(image, origin='lower', cmap='gray', extent = extent)
    ax[1, 0].plot(string_position['x0'], string_position['y0'], color='red')
    ax[1, 0].plot(string_position['x+'], string_position['y+'], 'r--')
    ax[1, 0].plot(string_position['x-'], string_position['y-'], 'r--')
    ax[1, 0].set_xlabel('Arcseconds')
    ax[1, 0].set_ylabel('Arcseconds')
    ax[1, 0].set_title('Cosmic String Position')

    # Plot the residuals
    ax[1, 1].imshow(np.abs(image - model), origin='lower', extent=extent, cmap='gray')
    #ax[1, 1].colorbar(label='Intensity')
    ax[1, 1].set_xlabel('Arcseconds')
    ax[1, 1].set_ylabel('Arcseconds')
    ax[1, 1].set_title('Residuals |data - model|')

    if saving:
        plt.imsave('fitted_lensed_sersic_galaxy.png', image, cmap='gray', origin='lower', extent=extent)

    plt.show()

def stringlensing_Model(params, model):

    pixel_scale = 0.05  # arcsec/pixel
    fwhm_arcsec = 0.8
    Ny, Nx = 256, 256
    model_square_ = model.reshape(Nx, Ny)
    Ie, re_arcsec, n, x0_pix, y0_pix, q, pos_angle_rad, P1, P2, pos_angle_string_rad, distane_center_string, RsRg = params
    tension = np.sqrt(P1**2 + P2**2) / (8 * np.pi)
    inclination = np.arctan(P2 / P1)

    sl.generateLensedSersicGalaxy_inplace(
        image = model_square_,
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

    model_square = sl.convolveWithGaussianPSF(
        model_square_,
        fwhm_arcsec=fwhm_arcsec,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    model = model_square.flatten()

    return model