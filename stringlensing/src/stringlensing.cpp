#define _USE_MATH_DEFINES
#define sec_in_rad 206264.80624709636

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h> 
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <functional>
#include <limits>
#include <iostream>

namespace py = pybind11;

// ================== ЯДРО: Sersic-профиль ==================
// при использовании matplotlib imshow позиционный угол будет отсчитываться в радианах против часовой стрелки от направления наверх.
// Предполагается, что симуляция производятся в коорданатах J2000 с полюсом мира наверху. Однако, в массиве это направление по убыванию координаты y.
double sersic_profile_point(
    double x,
    double y,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_rad,
    double pixel_scale_arcsec_per_pix
){
    double bn = 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n);
    double dx_pix = x - x0_pix;
    double dy_pix = y - y0_pix;

    double dx_pix_new = dy_pix * std::sin(pos_angle_rad) - dx_pix * std::cos(pos_angle_rad);
    double dy_pix_new = - dy_pix * std::cos(pos_angle_rad) - dx_pix * std::sin(pos_angle_rad);

    // расстояние в угловых секундах
    double R_arcsec = std::sqrt(dx_pix_new * dx_pix_new / (q * q) + dy_pix_new * dy_pix_new )
            * pixel_scale_arcsec_per_pix;

    double val = 0.0;
    if (re_arcsec > 0.0) {
        double ratio = R_arcsec / re_arcsec;
        double pow_term = std::pow(ratio, 1.0 / n);
        val = Ie * std::exp(-bn * (pow_term - 1.0));
    }
    return val;
}

void sersic_profile_inplace(
    double* img,
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_rad,
    double pixel_scale_arcsec_per_pix
) {
    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            img[y * Nx + x] = sersic_profile_point(
                static_cast<double>(x), static_cast<double>(y), 
                Ie,
                re_arcsec,
                n,
                x0_pix,
                y0_pix,
                q,
                pos_angle_rad,
                pixel_scale_arcsec_per_pix
            );
        }
    }
}

void lensed_sersic_profile_inplace(
    double* img_,
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_galaxy_rad,
    double tension,
    double inclination,
    double pos_angle_string_rad,
    double distane_center_string,
    double RsRg,
    double pixel_scale_arcsec_per_pix
) {
    double deficit = 8.0 * M_PI * tension * sec_in_rad / pixel_scale_arcsec_per_pix; // в пикселях
    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            // оригинальные координаты
            double dx_pix = static_cast<double>(x) - Nx / 2;
            double dy_pix = static_cast<double>(y) - Ny / 2;

            // Смещение задаёт положение струны относительно центра кадра.
            // Геометрию струны вычисляем в смещённых координатах, тогда как
            // профиль исходного объекта ниже берётся в координатах кадра.
            double x2 = dx_pix + distane_center_string * std::sin(pos_angle_string_rad);
            double y2 = dy_pix - distane_center_string * std::cos(pos_angle_string_rad);
            double dx_string = x2 * std::sin(pos_angle_string_rad) - y2 * std::cos(pos_angle_string_rad);
            double dy_string = x2 * std::cos(pos_angle_string_rad) + y2 * std::sin(pos_angle_string_rad);

            double RsY = RsRg / (1.0 + std::tan(inclination) * std::sin(dy_string * pixel_scale_arcsec_per_pix / sec_in_rad));
            double Re = deficit * (1.0 - RsY) * (std::cos(inclination) + std::sin(inclination) * dy_string * pixel_scale_arcsec_per_pix / sec_in_rad);

            double x_lensed_1 = x - Re * std::sin(pos_angle_string_rad) / 2;
            double y_lensed_1 = y + Re * std::cos(pos_angle_string_rad) / 2;
            double x_lensed_2 = x + Re * std::sin(pos_angle_string_rad) / 2;
            double y_lensed_2 = y - Re * std::cos(pos_angle_string_rad) / 2;

            if (dx_string >= -Re / 2){
                img_[y * Nx + x] += sersic_profile_point(
                    x_lensed_1, y_lensed_1,
                    Ie,
                    re_arcsec,
                    n,
                    x0_pix,
                    y0_pix,
                    q,
                    pos_angle_galaxy_rad,
                    pixel_scale_arcsec_per_pix
                );
            }
            if (dx_string <= Re / 2){
                img_[y * Nx + x] += sersic_profile_point(
                    x_lensed_2, y_lensed_2,
                    Ie,
                    re_arcsec,
                    n,
                    x0_pix,
                    y0_pix,
                    q,
                    pos_angle_galaxy_rad,
                    pixel_scale_arcsec_per_pix
                );
            }
        }
    }
}

void lensed_AB_sersic_profile_inplace(
    double* img_,
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_galaxy_rad,
    double ThetaE, // in arcseconds
    double dThetaEdXi, // arcseconds per radian of xi; converted to pixels below
    double pos_angle_string_rad,
    double distane_center_string,
    double pixel_scale_arcsec_per_pix
) {
    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            // оригинальные координаты
            double dx_pix = static_cast<double>(x) - Nx / 2;
            double dy_pix = static_cast<double>(y) - Ny / 2;

            // Смещение задаёт положение струны относительно центра кадра.
            // Геометрию струны вычисляем в смещённых координатах, тогда как
            // профиль исходного объекта ниже берётся в координатах кадра.
            double x2 = dx_pix + distane_center_string * std::sin(pos_angle_string_rad);
            double y2 = dy_pix - distane_center_string * std::cos(pos_angle_string_rad);
            double dx_string = x2 * std::sin(pos_angle_string_rad) - y2 * std::cos(pos_angle_string_rad);
            double dy_string = x2 * std::cos(pos_angle_string_rad) + y2 * std::sin(pos_angle_string_rad);

            double Re = ThetaE / pixel_scale_arcsec_per_pix + dThetaEdXi * dy_string / sec_in_rad;

            double x_lensed_1 = x - Re * std::sin(pos_angle_string_rad) / 2;
            double y_lensed_1 = y + Re * std::cos(pos_angle_string_rad) / 2;
            double x_lensed_2 = x + Re * std::sin(pos_angle_string_rad) / 2;
            double y_lensed_2 = y - Re * std::cos(pos_angle_string_rad) / 2;

            if (dx_string >= -Re / 2){
                img_[y * Nx + x] += sersic_profile_point(
                    x_lensed_1, y_lensed_1,
                    Ie,
                    re_arcsec,
                    n,
                    x0_pix,
                    y0_pix,
                    q,
                    pos_angle_galaxy_rad,
                    pixel_scale_arcsec_per_pix
                );
            }
            if (dx_string <= Re / 2){
                img_[y * Nx + x] += sersic_profile_point(
                    x_lensed_2, y_lensed_2,
                    Ie,
                    re_arcsec,
                    n,
                    x0_pix,
                    y0_pix,
                    q,
                    pos_angle_galaxy_rad,
                    pixel_scale_arcsec_per_pix
                );
            }
        }
    }
}

// ================== ЯДРО: Просто равномерный круг ==================

double uniform_profile_point(
    double x,
    double y,
    double Ie,
    double r0_arcsec,
    double x0_pix,
    double y0_pix,
    double pixel_scale_arcsec_per_pix
){
    double dx_pix = x - x0_pix;
    double dy_pix = y - y0_pix;

    // расстояние в угловых секундах
    double R_arcsec = std::sqrt(dx_pix * dx_pix + dy_pix * dy_pix) * pixel_scale_arcsec_per_pix;

    double val = 0.0;
    if (R_arcsec < r0_arcsec) {
        val = Ie;
    }
    return val;
}

void uniform_profile_inplace(
    double* img,
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double r0_arcsec,
    double x0_pix,
    double y0_pix,
    double pixel_scale_arcsec_per_pix
) {
    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            img[y * Nx + x] = uniform_profile_point(
                static_cast<double>(x), static_cast<double>(y), 
                Ie,
                r0_arcsec,
                x0_pix,
                y0_pix,
                pixel_scale_arcsec_per_pix
            );
        }
    }
}

void lensed_uniform_profile_inplace(
    double* img_,
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double r0_arcsec,
    double x0_pix,
    double y0_pix,
    double tension,
    double inclination,
    double pos_angle_string_rad,
    double distane_center_string,
    double RsRg,
    double pixel_scale_arcsec_per_pix
) {
    double deficit = 8.0 * M_PI * tension * sec_in_rad / pixel_scale_arcsec_per_pix; // в пикселях
    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            // оригинальные координаты
            double dx_pix = static_cast<double>(x) - Nx / 2;
            double dy_pix = static_cast<double>(y) - Ny / 2;

            // Смещение задаёт положение струны относительно центра кадра.
            // Геометрию струны вычисляем в смещённых координатах, тогда как
            // профиль исходного объекта ниже берётся в координатах кадра.
            double x2 = dx_pix + distane_center_string * std::sin(pos_angle_string_rad);
            double y2 = dy_pix - distane_center_string * std::cos(pos_angle_string_rad);
            double dx_string = x2 * std::sin(pos_angle_string_rad) - y2 * std::cos(pos_angle_string_rad);
            double dy_string = x2 * std::cos(pos_angle_string_rad) + y2 * std::sin(pos_angle_string_rad);

            double RsY = RsRg / (1.0 + std::tan(inclination) * std::sin(dy_string * pixel_scale_arcsec_per_pix / sec_in_rad));
            double Re = deficit * (1.0 - RsY) * (std::cos(inclination) + std::sin(inclination) * dy_string * pixel_scale_arcsec_per_pix / sec_in_rad);

            double x_lensed_1 = x - Re * std::sin(pos_angle_string_rad) / 2;
            double y_lensed_1 = y + Re * std::cos(pos_angle_string_rad) / 2;
            double x_lensed_2 = x + Re * std::sin(pos_angle_string_rad) / 2;
            double y_lensed_2 = y - Re * std::cos(pos_angle_string_rad) / 2;

            if (dx_string >= -Re / 2){
                img_[y * Nx + x] += uniform_profile_point(
                    x_lensed_1, y_lensed_1,
                    Ie,
                    r0_arcsec,
                    x0_pix,
                    y0_pix,
                    pixel_scale_arcsec_per_pix
                );
            }
            if (dx_string <= Re / 2){
                img_[y * Nx + x] += uniform_profile_point(
                    x_lensed_2, y_lensed_2,
                    Ie,
                    r0_arcsec,
                    x0_pix,
                    y0_pix,
                    pixel_scale_arcsec_per_pix
                );
            }
        }
    }
}

// ================== ЯДРО: свёртка с гауссом ==================

void gaussian_convolve_inplace(
    const double* img_in,
    double* img_out,
    ssize_t Ny,
    ssize_t Nx,
    double fwhm_arcsec,
    double pixel_scale_arcsec_per_pix
) {
    if (fwhm_arcsec <= 0.0) {
        // Нулевая/отрицательная FWHM — просто копируем
        std::copy(img_in, img_in + Ny * Nx, img_out);
        return;
    }

    // sigma в угл. секундах и в пикселях
    double sigma_arcsec = fwhm_arcsec / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    double sigma_pix = sigma_arcsec / pixel_scale_arcsec_per_pix;

    if (sigma_pix <= 0.0) {
        std::copy(img_in, img_in + Ny * Nx, img_out);
        return;
    }

    // Радиус ядра ~ 3 sigma
    int radius = static_cast<int>(std::ceil(3.0 * sigma_pix));
    int ksize = 2 * radius + 1;

    // 1D гауссово ядро
    std::vector<double> kernel(ksize);
    double norm = 0.0;
    for (int i = -radius; i <= radius; ++i) {
        double t = static_cast<double>(i) / sigma_pix;
        double val = std::exp(-0.5 * t * t);
        kernel[i + radius] = val;
        norm += val;
    }
    for (int i = 0; i < ksize; ++i) {
        kernel[i] /= norm;
    }

    // Временный буфер для промежуточной свёртки по X
    std::vector<double> tmp(Ny * Nx, 0.0);

    // Сначала свёртка по X
    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            double sum = 0.0;
            for (int dx = -radius; dx <= radius; ++dx) {
                ssize_t xx = x + dx;
                if (xx < 0 || xx >= Nx) continue;
                double w = kernel[dx + radius];
                sum += img_in[y * Nx + xx] * w;
            }
            tmp[y * Nx + x] = sum;
        }
    }

    // Затем свёртка по Y
    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            double sum = 0.0;
            for (int dy = -radius; dy <= radius; ++dy) {
                ssize_t yy = y + dy;
                if (yy < 0 || yy >= Ny) continue;
                double w = kernel[dy + radius];
                sum += tmp[yy * Nx + x] * w;
            }
            img_out[y * Nx + x] = sum;
        }
    }
}

// ================== Добавление шума и фона в картинки ==================

// Строго пуассоновский шум: каждый пиксель I -> k ~ Poisson( max(I, 0) )
void add_poisson_noise_inplace(py::array_t<double> img, 
                               double background,
                               std::uint64_t seed = 0)
{
    auto r = img.mutable_unchecked<2>();  // учитывает strides
    const ssize_t Ny = r.shape(0);
    const ssize_t Nx = r.shape(1);

    // Генератор случайных чисел
    std::mt19937_64 rng;
    if (seed == 0) {
        std::random_device rd;
        rng.seed(rd());
    } else {
        rng.seed(seed);
    }

    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            double I = r(y, x);

            // λ не может быть отрицательным
            double lambda = (I > 0.0) ? I : 0.0;

            // std::poisson_distribution ожидает параметр λ (double)
            std::poisson_distribution<long long> pois(lambda);

            long long k = pois(rng);  // количество "фотонов"

            // Записываем новое значение как целое число (в double)
            r(y, x) = static_cast<double>(k) + background;
        }
    }
}

// ================== Вычисление хи-квадрат ==================

double chi2_poisson(
    const py::array_t<double>& model,
    const py::array_t<double>& image,
    double min_error = 1e-3
) {
    auto m = model.unchecked<2>();  // быстрый доступ без копий
    auto d = image.unchecked<2>();

    if (m.shape(0) != d.shape(0) || m.shape(1) != d.shape(1)) {
        throw std::runtime_error("Model and image must have same shape.");
    }

    const ssize_t Ny = m.shape(0);
    const ssize_t Nx = m.shape(1);

    double chi2 = 0.0;

    for (ssize_t y = 0; y < Ny; ++y) {
        for (ssize_t x = 0; x < Nx; ++x) {
            double I = d(y, x);
            double M = m(y, x);

            // Пуассоновская ошибка: sigma = sqrt(I)
            // Защита от нуля/отрицательных значений
            double sigma = (I > 0.0) ? std::sqrt(I) : min_error;

            double diff = M - I;
            chi2 += (diff * diff) / (sigma * sigma);
        }
    }

    return chi2;
}

// ================= Минимизация ==================

// ---- Тип функции модели ----
using ModelFunc = std::function<void(const std::vector<double>& params,
                                     std::vector<double>& model)>;

struct LMResult {
    std::vector<double> params;
    double chi2;
    int iters;
    bool converged;

    std::vector<double> get_params() const { return params; }
    double get_chi2() const { return chi2; }
    int get_iters() const { return iters; }
    bool get_converged() const { return converged; }
};

// Обертка для Python функции
class PyModelFunc {
public:
    PyModelFunc(py::function func) : func_(func) {}
    
    void operator()(const std::vector<double>& params,
                    std::vector<double>& model) const {
        // Конвертируем vector в numpy array для Python
        py::array_t<double> py_params = py::array_t<double>(params.size(), params.data());
        
        // Создаем numpy array для результата
        py::array_t<double> py_model = py::array_t<double>(model.size());
        
        // Вызываем Python функцию
        py_model = func_(py_params, py_model);
        
        // Копируем результат обратно в C++ vector
        auto buf = py_model.request();
        double* ptr = static_cast<double*>(buf.ptr);
        std::copy(ptr, ptr + model.size(), model.begin());
    }
    
private:
    py::function func_;
};

// Простой решатель линейной системы A x = b методом Гаусса
// A: матрица N x N, разложенная по строкам, длина = N*N
// b: длина N
// x: длина N (результат)
bool solve_linear_system(std::vector<double>& A,
                         std::vector<double>& b,
                         std::vector<double>& x)
{
    const int N = static_cast<int>(b.size());
    x = b; // будем преобразовывать b в x

    // Прямой ход
    for (int k = 0; k < N; ++k) {
        // Поиск максимального по модулю элемента в столбце (частичный выбор главного элемента)
        double max_val = std::abs(A[k*N + k]);
        int max_row = k;
        for (int i = k + 1; i < N; ++i) {
            double val = std::abs(A[i*N + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        if (max_val < 1e-20) {
            return false; // Система вырожденная
        }

        // Меняем строки k и max_row в A и x
        if (max_row != k) {
            for (int j = 0; j < N; ++j) {
                std::swap(A[k*N + j], A[max_row*N + j]);
            }
            std::swap(x[k], x[max_row]);
        }

        // Нормализация и исключение
        double pivot = A[k*N + k];
        for (int j = k + 1; j < N; ++j) {
            A[k*N + j] /= pivot;
        }
        x[k] /= pivot;
        A[k*N + k] = 1.0;

        for (int i = k + 1; i < N; ++i) {
            double factor = A[i*N + k];
            if (std::abs(factor) < 1e-30) continue;

            for (int j = k + 1; j < N; ++j) {
                A[i*N + j] -= factor * A[k*N + j];
            }
            x[i] -= factor * x[k];
            A[i*N + k] = 0.0;
        }
    }

    // Обратный ход
    for (int i = N - 1; i >= 0; --i) {
        for (int j = i + 1; j < N; ++j) {
            x[i] -= A[i*N + j] * x[j];
        }
        // диагональ уже единичная
    }

    return true;
}

// ---- Собственно Levenberg–Marquardt для 2D изображений ----
LMResult levenberg_marquardt_2d(
    const ModelFunc& model_func,
    const std::vector<double>& image,
    int Ny, int Nx,
    const std::vector<double>& initial_params,
    int max_iters
) {
    const std::size_t Ndata = image.size();
    if (Ndata != static_cast<std::size_t>(Ny * Nx)) {
        throw std::runtime_error("image size != Ny * Nx");
    }

    const int Np = static_cast<int>(initial_params.size());
    if (Np <= 0) {
        throw std::runtime_error("number of parameters must be > 0");
    }

    // Текущие параметры
    std::vector<double> params = initial_params;

    // Буферы для модели, резидуалов и т.д.
    std::vector<double> model(Ndata);
    std::vector<double> model_pert(Ndata);
    std::vector<double> residuals(Ndata);      // r_i = (M_i - I_i)/sigma_i
    std::vector<double> residuals_pert(Ndata);

    // sigma_i (пуассоновские ошибки)
    std::vector<double> sigma(Ndata);
    for (std::size_t i = 0; i < Ndata; ++i) {
        double I = image[i];
        sigma[i] = (I > 0.0) ? std::sqrt(I) : 1.0;
    }

    // Начальное значение χ² и резидуалов
    model_func(params, model);
    double chi2 = 0.0;
    for (std::size_t i = 0; i < Ndata; ++i) {
        double r = (model[i] - image[i]) / sigma[i];
        residuals[i] = r;
        chi2 += r * r;
    }

    // Параметры алгоритма LM
    double lambda = 1e-2;      // начальный демпфирующий параметр
    const double lambda_up = 3.0;
    const double lambda_down = 0.3;

    const double eps_params = 1e-9;
    const double eps_chi2 = 1e-6;

    // Якобиан J: размер Ndata x Np
    std::vector<double> J(Ndata * Np);

    // Матрица A = J^T J и вектор g = J^T r
    std::vector<double> A(Np * Np);
    std::vector<double> g(Np);

    // Приращение параметров
    std::vector<double> delta(Np);

    bool converged = false;
    int iter = 0;

    for (; iter < max_iters; ++iter) {
        // --- Вычисляем якобиан по конечным разностям ---
        const double eps_fd = 1e-5;
        for (int k = 0; k < Np; ++k) {
            std::vector<double> params_pert = params;
            double step = eps_fd * (std::abs(params[k]) + 1.0);
            params_pert[k] += step;

            model_func(params_pert, model_pert);

            for (std::size_t i = 0; i < Ndata; ++i) {
                double r_plus = (model_pert[i] - image[i]) / sigma[i];
                double dr = (r_plus - residuals[i]) / step;  // d r_i / d p_k
                J[i * Np + k] = dr;
            }
        }

        // --- Строим A = J^T J и g = J^T r ---
        std::fill(A.begin(), A.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);

        for (std::size_t i = 0; i < Ndata; ++i) {
            const double ri = residuals[i];
            const double* Ji = &J[i * Np];
            for (int k = 0; k < Np; ++k) {
                g[k] += Ji[k] * ri;
                for (int l = 0; l < Np; ++l) {
                    A[k*Np + l] += Ji[k] * Ji[l];
                }
            }
        }

        // --- Добавляем демпфирование: (A + lambda * diag(A)) ---
        std::vector<double> A_damped = A;
        for (int k = 0; k < Np; ++k) {
            A_damped[k*Np + k] *= (1.0 + lambda);
        }

        // Решаем систему A_damped * delta = -g
        std::vector<double> rhs(Np);
        for (int k = 0; k < Np; ++k) {
            rhs[k] = -g[k];
        }

        if (!solve_linear_system(A_damped, rhs, delta)) {
            // Не удалось решить систему — выходим
            break;
        }

        // Проверка размера шага
        double max_rel_step = 0.0;
        for (int k = 0; k < Np; ++k) {
            double denom = std::abs(params[k]) + 1.0;
            double rel_step = std::abs(delta[k]) / denom;
            if (rel_step > max_rel_step) max_rel_step = rel_step;
        }
        if (max_rel_step < eps_params) {
            converged = true;
            break;
        }

        // Пробуем новые параметры
        std::vector<double> params_trial = params;
        for (int k = 0; k < Np; ++k) {
            params_trial[k] += delta[k];
        }

        model_func(params_trial, model_pert);

        double chi2_trial = 0.0;
        for (std::size_t i = 0; i < Ndata; ++i) {
            double r_trial = (model_pert[i] - image[i]) / sigma[i];
            chi2_trial += r_trial * r_trial;
        }

        // Принимаем/отклоняем шаг
        if (chi2_trial < chi2) {
            // Успех: принимаем
            params = std::move(params_trial);
            chi2 = chi2_trial;

            for (std::size_t i = 0; i < Ndata; ++i) {
                residuals[i] = (model_pert[i] - image[i]) / sigma[i];
                model[i] = model_pert[i];
            }

            // Уменьшаем λ, двигаемся к Гауссу–Ньютону
            lambda = std::max(lambda * lambda_down, 1e-12);

            // Проверка изменения χ²
            if (std::abs(chi2_trial - chi2) / (chi2 + 1e-12) < eps_chi2) {
                converged = true;

                break;
            }
        } else {
            // Неудачный шаг: увеличиваем λ (ближе к градиентному спуску)
            lambda = std::min(lambda * lambda_up, 1e12);
        }
    }

    LMResult res;
    res.params = params;
    res.chi2 = chi2;
    res.iters = iter;
    res.converged = converged;
    return res;
}

// ================= Генерация параметров =================

// ================== ОБЁРТКИ ДЛЯ PYBIND11 простого круга ==================

void generateUniformGalaxy_inplace(
    py::array_t<double> image,
    double Ie,
    double r0_arcsec,
    double x0_pix,
    double y0_pix,
    double pixel_scale_arcsec_per_pix
) {
    auto buf = image.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Image must be a 2D array");
    }
    ssize_t Ny = buf.shape[0];
    ssize_t Nx = buf.shape[1];

    double* ptr = static_cast<double*>(buf.ptr);

    uniform_profile_inplace(
        ptr, Ny, Nx,
        Ie, r0_arcsec,
        x0_pix, y0_pix,
        pixel_scale_arcsec_per_pix
    );
}

py::array_t<double> generateUniformGalaxy_new(
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double r0_arcsec,
    double x0_pix,
    double y0_pix,
    double pixel_scale_arcsec_per_pix
) {
    py::array_t<double> image({Ny, Nx});
    std::fill_n(static_cast<double*>(image.request().ptr), Ny * Nx, 0.0);
    auto buf = image.request();
    double* ptr = static_cast<double*>(buf.ptr);

    uniform_profile_inplace(
        ptr, Ny, Nx,
        Ie, r0_arcsec,
        x0_pix, y0_pix,
        pixel_scale_arcsec_per_pix
    );

    return image;
}

void generateLensedUniformGalaxy_inplace(
    py::array_t<double> image,
    double Ie,
    double r0_arcsec,
    double x0_pix,
    double y0_pix,
    double tension,
    double inclination,
    double pos_angle_string_rad,
    double distane_center_string,
    double RsRg,
    double pixel_scale_arcsec_per_pix
) {
    auto buf = image.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Image must be a 2D array");
    }
    ssize_t Ny = buf.shape[0];
    ssize_t Nx = buf.shape[1];

    double* ptr = static_cast<double*>(buf.ptr);

    lensed_uniform_profile_inplace(
        ptr, Ny, Nx,
        Ie, r0_arcsec,
        x0_pix, y0_pix,
        tension, inclination,
        pos_angle_string_rad,
        distane_center_string,
        RsRg,
        pixel_scale_arcsec_per_pix
    );
}

py::array_t<double> generateLensedUniformGalaxy_new(
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double r0_arcsec,
    double x0_pix,
    double y0_pix,
    double tension,
    double inclination,
    double pos_angle_string_rad,
    double distane_center_string,
    double RsRg,
    double pixel_scale_arcsec_per_pix
) {
    py::array_t<double> image({Ny, Nx});
    std::fill_n(static_cast<double*>(image.request().ptr), Ny * Nx, 0.0);
    auto buf = image.request();
    double* ptr = static_cast<double*>(buf.ptr);

    lensed_uniform_profile_inplace(
        ptr, Ny, Nx,
        Ie, r0_arcsec,
        x0_pix, y0_pix,
        tension, inclination,
        pos_angle_string_rad,
        distane_center_string,
        RsRg,
        pixel_scale_arcsec_per_pix
    );

    return image;
}

// ================== ОБЁРТКИ ДЛЯ PYBIND11 серсика ==================

// 1) In-place генерация Sersic-профиля в переданном массиве
void generateSersicGalaxy_inplace(
    py::array_t<double> image,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_rad,
    double pixel_scale_arcsec_per_pix
) {
    auto buf = image.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Image must be a 2D array");
    }
    ssize_t Ny = buf.shape[0];
    ssize_t Nx = buf.shape[1];

    double* ptr = static_cast<double*>(buf.ptr);

    sersic_profile_inplace(
        ptr, Ny, Nx,
        Ie, re_arcsec, n,
        x0_pix, y0_pix, q, pos_angle_rad,
        pixel_scale_arcsec_per_pix
    );
}

// 2) Вариант, который СОЗДАЁТ новый массив и возвращает его
py::array_t<double> generateSersicGalaxy_new(
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_rad,
    double pixel_scale_arcsec_per_pix
) {
    py::array_t<double> image({Ny, Nx});
    std::fill_n(static_cast<double*>(image.request().ptr), Ny * Nx, 0.0);
    auto buf = image.request();
    double* ptr = static_cast<double*>(buf.ptr);

    sersic_profile_inplace(
        ptr, Ny, Nx,
        Ie, re_arcsec, n,
        x0_pix, y0_pix,q, pos_angle_rad,
        pixel_scale_arcsec_per_pix
    );

    return image;
}

// 1) In-place генерация линзированного Sersic-профиля в переданном массиве
void generateLensedSersicGalaxy_inplace(
    py::array_t<double> image,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_galaxy_rad,
    double tension,
    double inclination,
    double pos_angle_string_rad,
    double distane_center_string,
    double RsRg,
    double pixel_scale_arcsec_per_pix
) {
    auto buf = image.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Image must be a 2D array");
    }
    ssize_t Ny = buf.shape[0];
    ssize_t Nx = buf.shape[1];

    double* ptr = static_cast<double*>(buf.ptr);

    lensed_sersic_profile_inplace(
        ptr, Ny, Nx,
        Ie, re_arcsec, n,
        x0_pix, y0_pix, q, pos_angle_galaxy_rad,
        tension, inclination,
        pos_angle_string_rad,
        distane_center_string,
        RsRg,
        pixel_scale_arcsec_per_pix
    );
}

// 2) Вариант, который СОЗДАЁТ новый массив и возвращает его
py::array_t<double> generateLensedABSersicGalaxy_new(
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_galaxy_rad,
    double ThetaE, // in arcseconds
    double dThetaEdXi, // arcseconds per radian of xi; converted to pixels below
    double pos_angle_string_rad,
    double distane_center_string,
    double pixel_scale_arcsec_per_pix
) {
    py::array_t<double> image({Ny, Nx});
    std::fill_n(static_cast<double*>(image.request().ptr), Ny * Nx, 0.0);
    auto buf = image.request();
    double* ptr = static_cast<double*>(buf.ptr);

    lensed_AB_sersic_profile_inplace(
        ptr, Ny, Nx,
        Ie, re_arcsec, n,
        x0_pix, y0_pix, q, pos_angle_galaxy_rad,
        ThetaE, dThetaEdXi,
        pos_angle_string_rad,
        distane_center_string,
        pixel_scale_arcsec_per_pix
    );

    return image;
}

py::array_t<double> generateLensedSersicGalaxy_new(
    ssize_t Ny,
    ssize_t Nx,
    double Ie,
    double re_arcsec,
    double n,
    double x0_pix,
    double y0_pix,
    double q,
    double pos_angle_galaxy_rad,
    double tension,
    double inclination,
    double pos_angle_string_rad,
    double distane_center_string,
    double RsRg,
    double pixel_scale_arcsec_per_pix
) {
    py::array_t<double> image({Ny, Nx});
    std::fill_n(static_cast<double*>(image.request().ptr), Ny * Nx, 0.0);
    auto buf = image.request();
    double* ptr = static_cast<double*>(buf.ptr);

    lensed_sersic_profile_inplace(
        ptr, Ny, Nx,
        Ie, re_arcsec, n,
        x0_pix, y0_pix, q, pos_angle_galaxy_rad,
        tension, inclination,
        pos_angle_string_rad,
        distane_center_string,
        RsRg,
        pixel_scale_arcsec_per_pix
    );

    return image;
}

// 3) Свёртка любой карты яркости с гауссовой PSF (FWHM в arcsec)
py::array_t<double> convolveWithGaussianPSF(
    py::array_t<double> image,
    double fwhm_arcsec,
    double pixel_scale_arcsec_per_pix
) {
    auto buf = image.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Image must be a 2D array");
    }
    ssize_t Ny = buf.shape[0];
    ssize_t Nx = buf.shape[1];

    // На всякий случай копируем вход в отдельный буфер (чтобы не портить оригинал)
    std::vector<double> img_in(Ny * Nx);
    {
        double* src = static_cast<double*>(buf.ptr);
        std::copy(src, src + Ny * Nx, img_in.begin());
    }

    // Выходной NumPy-массив
    py::array_t<double> out({Ny, Nx});
    auto buf_out = out.request();
    double* img_out = static_cast<double*>(buf_out.ptr);

    gaussian_convolve_inplace(
        img_in.data(),
        img_out,
        Ny,
        Nx,
        fwhm_arcsec,
        pixel_scale_arcsec_per_pix
    );

    return out;
}

// 4) Обертка фитирования для вызова из Python
LMResult levenberg_marquardt_2d_py(
    py::function py_model_func,
    py::array_t<double> image_array,
    int Ny, int Nx,
    py::array_t<double> initial_params_array,
    int max_iters
) {
    // Конвертируем numpy arrays в C++ vectors
    auto image_buf = image_array.request();
    auto params_buf = initial_params_array.request();
    
    if (image_buf.ndim != 1) {
        throw std::runtime_error("Image must be 1D array (flattened)");
    }
    
    if (params_buf.ndim != 1) {
        throw std::runtime_error("Initial params must be 1D array");
    }
    
    std::vector<double> image(
        static_cast<double*>(image_buf.ptr),
        static_cast<double*>(image_buf.ptr) + image_buf.size
    );
    
    std::vector<double> initial_params(
        static_cast<double*>(params_buf.ptr),
        static_cast<double*>(params_buf.ptr) + params_buf.size
    );
    
    // Создаем обертку для Python функции
    PyModelFunc cpp_model_func(py_model_func);
    
    // Вызываем C++ алгоритм
    return levenberg_marquardt_2d(
        cpp_model_func,
        image,
        Ny, Nx,
        initial_params,
        max_iters
    );
}

// 4) Обертка фитирования для вызова из Python с поддержкой 2D numpy arrays
LMResult levenberg_marquardt_2d_py_2d(
    py::function py_model_func,
    py::array_t<double> image_array_2d,
    py::array_t<double> initial_params_array,
    int max_iters
) {
    auto image_buf = image_array_2d.request();
    auto params_buf = initial_params_array.request();
    
    if (image_buf.ndim != 2) {
        throw std::runtime_error("Image must be 2D array");
    }
    
    // Получаем размеры
    int Ny = image_buf.shape[0];
    int Nx = image_buf.shape[1];
    
    // Преобразуем 2D массив в 1D vector
    double* image_ptr = static_cast<double*>(image_buf.ptr);
    std::vector<double> image(image_ptr, image_ptr + Ny * Nx);
    
    std::vector<double> initial_params(
        static_cast<double*>(params_buf.ptr),
        static_cast<double*>(params_buf.ptr) + params_buf.size
    );
    
    PyModelFunc cpp_model_func(py_model_func);
    
    return levenberg_marquardt_2d(
        cpp_model_func,
        image,
        Ny, Nx,
        initial_params,
        max_iters
    );
}

// Дальнейшие планы:
// 3) Исправить алгоритм минимизации Levenberg-Marquardt method
// 6) Добавить генерацию параметров для тестирования точности алгоритма
// 7) Распараллелить вычисления с помощью OpenMP

// ================== МОДУЛЬ PYBIND11 ==================

PYBIND11_MODULE(stringlensing, sl) {
    sl.doc() = "Operational module for fitting later in lmfit";

    // Константа
    sl.attr("version") = "1.1";

    sl.def(
        "generateUniformGalaxy_inplace",
        &generateUniformGalaxy_inplace,
        "Fill a given 2D NumPy array with a uniform galaxy profile",
        py::arg("image"),
        py::arg("Ie"),
        py::arg("r0_arcsec"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateUniformGalaxy_new",
        &generateUniformGalaxy_new,
        "Create and return a new uniform galaxy image",
        py::arg("Ny"),
        py::arg("Nx"),
        py::arg("Ie"),
        py::arg("r0_arcsec"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateLensedUniformGalaxy_inplace",
        &generateLensedUniformGalaxy_inplace,
        "Fill a given 2D NumPy array with a lensed uniform galaxy profile",
        py::arg("image"),
        py::arg("Ie"),
        py::arg("r0_arcsec"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("tension"),
        py::arg("inclination"),
        py::arg("pos_angle_string_rad"),
        py::arg("distane_center_string"),
        py::arg("RsRg"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateLensedUniformGalaxy_new",
        &generateLensedUniformGalaxy_new,
        "Create and return a new lensed uniform galaxy image",
        py::arg("Ny"),
        py::arg("Nx"),
        py::arg("Ie"),
        py::arg("r0_arcsec"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("tension"),
        py::arg("inclination"),
        py::arg("pos_angle_string_rad"),
        py::arg("distane_center_string"),
        py::arg("RsRg"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateSersicGalaxy_inplace",
        &generateSersicGalaxy_inplace,
        "Fill a given 2D NumPy array with a Sersic galaxy profile",
        py::arg("image"),
        py::arg("Ie"),
        py::arg("re_arcsec"),
        py::arg("n"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("q"),
        py::arg("pos_angle_rad"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateSersicGalaxy_new",
        &generateSersicGalaxy_new,
        "Create and return a new Sersic galaxy image",
        py::arg("Ny"),
        py::arg("Nx"),
        py::arg("Ie"),
        py::arg("re_arcsec"),
        py::arg("n"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("q"),
        py::arg("pos_angle_rad"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateLensedSersicGalaxy_inplace",
        &generateLensedSersicGalaxy_inplace,
        "Fill a given 2D NumPy array with a lensed Sersic galaxy profile",
        py::arg("image"),
        py::arg("Ie"),
        py::arg("re_arcsec"),
        py::arg("n"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("q"),
        py::arg("pos_angle_galaxy_rad"),
        py::arg("tension"),
        py::arg("inclination"),
        py::arg("pos_angle_string_rad"),
        py::arg("distane_center_string"),
        py::arg("RsRg"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateLensedABSersicGalaxy_new",
        &generateLensedABSersicGalaxy_new,
        "Create and return a new lensed AB Sersic galaxy image",
        py::arg("Ny"),
        py::arg("Nx"),
        py::arg("Ie"),
        py::arg("re_arcsec"),
        py::arg("n"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("q"),
        py::arg("pos_angle_galaxy_rad"),
        py::arg("ThetaE"),
        py::arg("dThetaEdXi"),
        py::arg("pos_angle_string_rad"),
        py::arg("distane_center_string"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "generateLensedSersicGalaxy_new",
        &generateLensedSersicGalaxy_new,
        "Create and return a new lensed Sersic galaxy image",
        py::arg("Ny"),
        py::arg("Nx"),
        py::arg("Ie"),
        py::arg("re_arcsec"),
        py::arg("n"),
        py::arg("x0_pix"),
        py::arg("y0_pix"),
        py::arg("q"),
        py::arg("pos_angle_galaxy_rad"),
        py::arg("tension"),
        py::arg("inclination"),
        py::arg("pos_angle_string_rad"),
        py::arg("distane_center_string"),
        py::arg("RsRg"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "convolveWithGaussianPSF",
        &convolveWithGaussianPSF,
        "Convolve a 2D image with a Gaussian PSF (given FWHM in arcsec)",
        py::arg("image"),
        py::arg("fwhm_arcsec"),
        py::arg("pixel_scale_arcsec_per_pix")
    );

    sl.def(
        "add_poisson_noise_inplace",
        &add_poisson_noise_inplace,
        py::arg("img"),
        py::arg("background") = 0.0,
        py::arg("seed") = 0,
        "Add strictly Poisson noise to image in-place"
    );

    sl.def("chi2_poisson", &chi2_poisson,
          py::arg("model"), 
          py::arg("image"), 
          py::arg("min_error"),
          "Compute chi^2 with Poisson errors sqrt(image)");

    // Экспортируем структуру LMResult
    py::class_<LMResult>(sl, "LMResult")
        .def(py::init<>())
        .def_readwrite("params", &LMResult::params)
        .def_readwrite("chi2", &LMResult::chi2)
        .def_readwrite("iters", &LMResult::iters)
        .def_readwrite("converged", &LMResult::converged)
        .def("get_params", &LMResult::get_params)
        .def("get_chi2", &LMResult::get_chi2)
        .def("get_iters", &LMResult::get_iters)
        .def("get_converged", &LMResult::get_converged)
        .def("__repr__", [](const LMResult& r) {
            std::string repr = "LMResult(\n";
            repr += "  params=[";
            for (size_t i = 0; i < r.params.size(); ++i) {
                repr += std::to_string(r.params[i]);
                if (i < r.params.size() - 1) repr += ", ";
            }
            repr += "],\n";
            repr += "  chi2=" + std::to_string(r.chi2) + ",\n";
            repr += "  iters=" + std::to_string(r.iters) + ",\n";
            repr += "  converged=" + std::string(r.converged ? "True" : "False") + "\n";
            repr += ")";
            return repr;
        });
    
    // Экспортируем функцию с разными сигнатурами
    sl.def("levenberg_marquardt_2d", 
          &levenberg_marquardt_2d_py,
          "Levenberg-Marquardt optimization for 2D images",
          py::arg("model_func"),
          py::arg("image"),           // 1D flattened array
          py::arg("Ny"),              // height (rows)
          py::arg("Nx"),              // width (columns)
          py::arg("initial_params"),  // 1D array
          py::arg("max_iters") = 100);
    
    sl.def("levenberg_marquardt_2d_py_2d", 
          &levenberg_marquardt_2d_py_2d,
          "Levenberg-Marquardt optimization for 2D images (2D array input)",
          py::arg("model_func"),
          py::arg("image"),           // 2D array
          py::arg("initial_params"),  // 1D array
          py::arg("max_iters") = 100);
}