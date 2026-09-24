from __future__ import annotations

#multiprocessing version 1
'''
from concurrent.futures import ProcessPoolExecutor
'''

#multiprocessing version 2 
import ray

import os
import subprocess
import argparse
import json
import math
import sys
import traceback
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib-fit-quality"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Parameters, minimize
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr, spearmanr

import stringlensing as sl
import stringlensing_utilities as sl_utils



# =============================================================================
# Конфигурация и абстракции
# =============================================================================


Array2D = np.ndarray
ParamDict = Dict[str, float]
EnvDict = Dict[str, Any]


@dataclass
class ParameterSpec:
    name: str
    min_value: float
    max_value: float
    truth_sampler: Callable[[np.random.Generator, EnvDict], float]
    init_sampler: Optional[Callable[[np.random.Generator, EnvDict], float]] = None
    vary: bool = True

    def sample_truth(self, rng: np.random.Generator, env: EnvDict) -> float:
        value = float(self.truth_sampler(rng, env))
        return float(np.clip(value, self.min_value, self.max_value))

    def sample_init(self, rng: np.random.Generator, env: EnvDict) -> float:
        if self.init_sampler is None:
            value = rng.uniform(self.min_value, self.max_value)
        else:
            value = float(self.init_sampler(rng, env))
        return float(np.clip(value, self.min_value, self.max_value))


@dataclass
class FitSettings:
    methods_chain: List[str] = field(default_factory=lambda: ["least_squares", "nelder", "leastsq"])
    max_retries: int = 3
    nan_policy: str = "raise"
    calc_covar: bool = True
    max_nfev: Optional[int] = 5000
    residual_sigma_source: str = "truth_noise"  # truth_noise | unity
    require_positive_stderr: bool = True


@dataclass
class PipelineConfig:
    output_dir: str = "fit_quality_project"
    n_samples: int = 50
    random_seed: int = 42
    example_plots_to_save: int = 6
    save_per_sample_fit_json: bool = True
    save_example_images: bool = True
    image_dtype: str = "float32"


@dataclass
class RunSummary:
    n_total: int
    n_success: int
    n_failed: int
    failed_sample_ids: List[str]


class ModelAdapter:
    """Обертка для внешней функции вида:

    external_func(nx, ny, *model_params, pixel_scale, noise_sigma, psf_sigma, ...)

    Для вашего реального случая достаточно заменить `external_func` и списки имен
    параметров. Весь остальной pipeline останется тем же.
    """

    def __init__(
        self,
        external_func: Callable[..., Array2D],
        model_param_names: Sequence[str],
        env_param_names: Sequence[str],
    ) -> None:
        self.external_func = external_func
        self.model_param_names = list(model_param_names)
        self.env_param_names = list(env_param_names)

    def __call__(self, model_params: ParamDict, env: EnvDict, rng: Optional[np.random.Generator] = None) -> Array2D:
        nx = int(env["nx"])
        ny = int(env["ny"])
        ordered_model_params = [model_params[name] for name in self.model_param_names]
        kwargs = {name: env[name] for name in self.env_param_names if name not in {"nx", "ny"}}
        if "rng" in self.external_func.__code__.co_varnames:
            kwargs["rng"] = rng
        image = self.external_func(nx, ny, *ordered_model_params, **kwargs)
        return np.asarray(image, dtype=np.float64)


# =============================================================================
# Пример модели: сумма двух 2D-гауссиан + фон = 11 параметров модели
# =============================================================================


def elliptical_gaussian_2d(
    x_arcsec: Array2D,
    y_arcsec: Array2D,
    amp: float,
    x0: float,
    y0: float,
    sx: float,
    sy: float,
) -> Array2D:
    sx = max(float(sx), 1e-6)
    sy = max(float(sy), 1e-6)
    return amp * np.exp(-0.5 * (((x_arcsec - x0) / sx) ** 2 + ((y_arcsec - y0) / sy) ** 2))


def toy_external_two_gaussians(
    nx: int,
    ny: int,
    amp1: float,
    x01: float,
    y01: float,
    sx1: float,
    sy1: float,
    amp2: float,
    x02: float,
    y02: float,
    sx2: float,
    sy2: float,
    background: float,
    pixel_scale: float,
    noise_sigma: float,
    psf_sigma: float,
    rng: Optional[np.random.Generator] = None,
) -> Array2D:
    """Тестовая внешняя функция с той же логикой интерфейса, что и у пользователя.

    Параметры модели (11):
      amp1, x01, y01, sx1, sy1,
      amp2, x02, y02, sx2, sy2,
      background

    Параметры окружения (5):
      nx, ny, pixel_scale, noise_sigma, psf_sigma
    """
    x = (np.arange(nx) - (nx - 1) / 2.0) * pixel_scale
    y = (np.arange(ny) - (ny - 1) / 2.0) * pixel_scale
    xx, yy = np.meshgrid(x, y, indexing="ij")

    image = (
        elliptical_gaussian_2d(xx, yy, amp1, x01, y01, sx1, sy1)
        + elliptical_gaussian_2d(xx, yy, amp2, x02, y02, sx2, sy2)
        + background
    )

    psf_sigma_pix = max(float(psf_sigma) / float(pixel_scale), 0.0)
    if psf_sigma_pix > 0:
        image = gaussian_filter(image, sigma=psf_sigma_pix, mode="nearest")

    if noise_sigma > 0:
        if rng is None:
            rng = np.random.default_rng()
        image = image + rng.normal(loc=0.0, scale=noise_sigma, size=image.shape)

    return image

# =============================================================================
# Реальная модель: галактика с профилем серсика + космическая струна = 11 параметров модели
# =============================================================================

def real_external_galaxy_plus_string(
    ny: int,
    nx: int,
    Ie: float,
    re_arcsec: float,
    n: float,
    x0_pix: float,
    y0_pix: float,
    q: float,
    pos_angle_rad: float,
    ThetaE: float,
    dThetaEdXi: float,
    pos_angle_string_rad: float,
    distane_center_string: float,
    pixel_scale: float,
    noise_sigma: float,
    psf_sigma: float,
    rng: Optional[np.random.Generator] = None,
) -> Array2D:

    '''
    ThetaE = 2.0
    dThetaEdXi = 0.0
    '''

    img_11 = sl.generateLensedABSersicGalaxy_new(
        ny, nx,
        Ie=Ie,
        re_arcsec=re_arcsec,
        n=n,
        x0_pix=x0_pix + nx / 2.0,
        y0_pix=y0_pix + ny / 2.0,
        q=q,
        pos_angle_galaxy_rad = pos_angle_rad,
        ThetaE = ThetaE,
        dThetaEdXi = dThetaEdXi,
        pos_angle_string_rad = pos_angle_string_rad,
        distane_center_string = distane_center_string,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    if noise_sigma > 0:
        if rng is None:
            rng = np.random.default_rng()
        seed = abs(int(rng.normal(loc=0.0, scale=1000, size=1)[0]))
        sl.add_poisson_noise_inplace(
            img_11, 
            background=0.0, 
            seed=seed
        )

    img_12 = sl.convolveWithGaussianPSF(
        img_11,
        fwhm_arcsec=psf_sigma,
        pixel_scale_arcsec_per_pix=pixel_scale
    )

    return img_12

# =============================================================================
# Спецификация тетсовой задачи (11 параметров)
# =============================================================================

def build_toy_parameter_specs(env: EnvDict) -> List[ParameterSpec]:
    fov_x = env["nx"] * env["pixel_scale"]
    fov_y = env["ny"] * env["pixel_scale"]
    x_lim = 0.35 * fov_x
    y_lim = 0.35 * fov_y

    def uniform_sampler(lo: float, hi: float) -> Callable[[np.random.Generator, EnvDict], float]:
        return lambda rng, _env: rng.uniform(lo, hi)

    return [
        ParameterSpec("amp1", 10.0, 80.0, uniform_sampler(20.0, 70.0)),
        ParameterSpec("x01", -x_lim, x_lim, uniform_sampler(-0.25 * fov_x, 0.25 * fov_x)),
        ParameterSpec("y01", -y_lim, y_lim, uniform_sampler(-0.25 * fov_y, 0.25 * fov_y)),
        ParameterSpec("sx1", 0.25, 2.50, uniform_sampler(0.35, 1.50)),
        ParameterSpec("sy1", 0.25, 2.50, uniform_sampler(0.35, 1.50)),
        ParameterSpec("amp2", 10.0, 80.0, uniform_sampler(20.0, 70.0)),
        ParameterSpec("x02", -x_lim, x_lim, uniform_sampler(-0.25 * fov_x, 0.25 * fov_x)),
        ParameterSpec("y02", -y_lim, y_lim, uniform_sampler(-0.25 * fov_y, 0.25 * fov_y)),
        ParameterSpec("sx2", 0.25, 2.50, uniform_sampler(0.35, 1.50)),
        ParameterSpec("sy2", 0.25, 2.50, uniform_sampler(0.35, 1.50)),
        ParameterSpec("background", 0.0, 10.0, uniform_sampler(0.2, 3.0)),
    ]

# =============================================================================
# Спецификация реальной задачи (11 параметров)
# =============================================================================

def build_real_parameter_specs(env: EnvDict) -> List[ParameterSpec]:
    fov_x = env["nx"] / 2 #* env["pixel_scale"]
    fov_y = env["ny"] / 2 #* env["pixel_scale"]
    x_lim = 0.8 * fov_x
    y_lim = 0.8 * fov_y

    def uniform_sampler(lo: float, hi: float) -> Callable[[np.random.Generator, EnvDict], float]:
        return lambda rng, _env: rng.uniform(lo, hi)

    return [
        ParameterSpec("Ie", 0.5, 60.0, uniform_sampler(1.0, 50.0)),
        ParameterSpec("re_arcsec", 1.0, 7.0, uniform_sampler(2.0, 6.0)),
        ParameterSpec("n", 0.5, 4.5, uniform_sampler(1.0, 4.0)),
        ParameterSpec("x0_pix", -x_lim, x_lim, uniform_sampler(-0.5 * fov_x, 0.5 * fov_x)),
        ParameterSpec("y0_pix", -y_lim, y_lim, uniform_sampler(-0.5 * fov_y, 0.5 * fov_y)),
        ParameterSpec("q", 0.1, 1.0, uniform_sampler(0.2, 0.9)),
        ParameterSpec("pos_angle_rad", 0.0, np.pi, uniform_sampler(0, np.pi)),
        ParameterSpec("ThetaE", 0.1, 5.0, uniform_sampler(0.5, 3.5)),
        ParameterSpec("dThetaEdXi", 0.0, 206265, uniform_sampler(0.0, 206265)),
        ParameterSpec("pos_angle_string_rad", 0.0, np.pi, uniform_sampler(0, np.pi)),
        ParameterSpec("distane_center_string", 0.0, 3.0, uniform_sampler(0.0, 3.0)),
    ]

# =============================================================================
# Общие утилиты
# =============================================================================


def ensure_dirs(base_dir: Path) -> Dict[str, Path]:
    dirs = {
        "root": base_dir,
        "synthetic_images": base_dir / "synthetic" / "images",
        "synthetic_truth": base_dir / "synthetic" / "truth",
        "fit_json": base_dir / "fit_results" / "json",
        "plots": base_dir / "analysis" / "plots",
        "tables": base_dir / "analysis" / "tables",
        "examples": base_dir / "analysis" / "examples",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs



def to_jsonable(obj: Any) -> Any:
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    return obj


# =============================================================================
# Блок 1. Генерация синтетики
# =============================================================================


def sample_truth_parameters(
    specs: Sequence[ParameterSpec],
    env: EnvDict,
    rng: np.random.Generator,
) -> ParamDict:
    values = {spec.name: spec.sample_truth(rng, env) for spec in specs}
    return values



def maybe_relabel_two_sources(params: ParamDict) -> ParamDict:
    """Устраняет перестановочную неоднозначность двух гауссиан.

    Сортируем компоненты по x-координате. Это важно для корректного сравнения
    true vs fitted, чтобы фит не считался плохим только из-за swap двух гауссиан.
    """
    if {"amp1", "x01", "y01", "sx1", "sy1", "amp2", "x02", "y02", "sx2", "sy2"}.issubset(params.keys()):
        g1 = {k[:-1]: params[k] for k in ["amp1", "x01", "y01", "sx1", "sy1"]}
        g2 = {k[:-1]: params[k] for k in ["amp2", "x02", "y02", "sx2", "sy2"]}
        if g2["x0"] < g1["x0"]:
            out = dict(params)
            for base in ["amp", "x0", "y0", "sx", "sy"]:
                out[f"{base}1"] = g2[base]
                out[f"{base}2"] = g1[base]
            return out
    return params



def generate_synthetic_dataset(
    adapter: ModelAdapter,
    specs: Sequence[ParameterSpec],
    fixed_env: EnvDict,
    pipeline_cfg: PipelineConfig,
) -> pd.DataFrame:
    base_dir = Path(pipeline_cfg.output_dir)
    dirs = ensure_dirs(base_dir)
    rng = np.random.default_rng(pipeline_cfg.random_seed)

    manifest_rows: List[Dict[str, Any]] = []

    for idx in range(pipeline_cfg.n_samples):
        sample_id = f"sample_{idx:05d}"
        true_params = maybe_relabel_two_sources(sample_truth_parameters(specs, fixed_env, rng))
        image = adapter(true_params, fixed_env, rng=rng).astype(pipeline_cfg.image_dtype)

        image_path = dirs["synthetic_images"] / f"{sample_id}.npy"
        truth_path = dirs["synthetic_truth"] / f"{sample_id}.json"
        np.save(image_path, image)

        truth_payload = {
            "sample_id": sample_id,
            "true_params": true_params,
            "env": fixed_env,
            "image_path": str(image_path),
        }
        truth_path.write_text(json.dumps(to_jsonable(truth_payload), indent=2), encoding="utf-8")

        row: Dict[str, Any] = {
            "sample_id": sample_id,
            "image_path": str(image_path),
            "truth_json_path": str(truth_path),
        }
        row.update({f"true_{k}": v for k, v in true_params.items()})
        row.update({f"env_{k}": v for k, v in fixed_env.items()})
        manifest_rows.append(row)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(base_dir / "synthetic_manifest.csv", index=False)
    return manifest


# =============================================================================
# Блок 2. Фитирование
# =============================================================================


def build_lmfit_parameters(
    specs: Sequence[ParameterSpec],
    env: EnvDict,
    rng: np.random.Generator,
) -> Parameters:
    params = Parameters()
    for spec in specs:
        params.add(
            name=spec.name,
            value=spec.sample_init(rng, env),
            vary=spec.vary,
            min=spec.min_value,
            max=spec.max_value,
        )
    return params



def lmfit_params_to_dict(params: Parameters) -> ParamDict:
    return {name: float(par.value) for name, par in params.items()}



def residual_function(
    lm_params: Parameters,
    image_data: Array2D,
    adapter: ModelAdapter,
    render_env: EnvDict,
    sigma_for_weights: Array2D,
) -> np.ndarray:
    model_params = lmfit_params_to_dict(lm_params)
    model_image = adapter(model_params, render_env, rng=None)
    resid = model_image - image_data
    resid = resid / sigma_for_weights
    return np.asarray(resid.ravel(), dtype=np.float64)



def fit_result_has_valid_uncertainties(result: Any, params: Parameters) -> Tuple[bool, str]:
    if getattr(result, "covar", None) is None:
        return False, "covariance matrix is None"
    if not bool(getattr(result, "errorbars", False)):
        return False, "errorbars flag is False"
    for name, par in params.items():
        if not getattr(par, "vary", True):
            continue
        stderr = getattr(par, "stderr", None)
        if stderr is None:
            return False, f"stderr is None for parameter '{name}'"
        if not np.isfinite(stderr):
            return False, f"stderr is not finite for parameter '{name}'"
        if stderr <= 0:
            return False, f"stderr <= 0 for parameter '{name}'"
    return True, "ok"



def run_fit_once(
    image_data: Array2D,
    adapter: ModelAdapter,
    specs: Sequence[ParameterSpec],
    sample_env: EnvDict,
    fit_cfg: FitSettings,
    rng: np.random.Generator,
) -> Any:
    params = build_lmfit_parameters(specs, sample_env, rng)
    render_env = dict(sample_env)
    render_env["noise_sigma"] = 0.0
    
    # gaussian noise
    sigma_value = float(sample_env.get("noise_sigma", 0.0)) if fit_cfg.residual_sigma_source == "truth_noise" else 1.0
    sigma = np.full_like(image_data, fill_value=sigma_value, dtype=np.float64)

    # poisson noise
    sigma = np.sqrt(image_data.clip(min=0.0) + sigma**2)

    result = None
    for i, method in enumerate(fit_cfg.methods_chain):
        kwargs = dict(
            kws={
                "image_data": image_data,
                "adapter": adapter,
                "render_env": render_env,
                "sigma_for_weights": sigma,
            },
            method=method,
            nan_policy=fit_cfg.nan_policy,
            calc_covar=fit_cfg.calc_covar,
            max_nfev=fit_cfg.max_nfev,
        )
        result = minimize(residual_function, params, **kwargs)
        params = result.params.copy()
    return result


@ray.remote
def run_fit_retrying(rec, adapter, specs, fit_cfg, synthetic_manifest, rng) -> Tuple[Optional[Exception], Dict[str, Any]]:
    sample_id = rec["sample_id"]
    image_data = np.load(rec["image_path"]).astype(np.float64)
    sample_env = {
        col.replace("env_", ""): rec[col]
        for col in synthetic_manifest.columns
        if col.startswith("env_")
    }

    best_payload: Optional[Dict[str, Any]] = None
    last_exception = None
    
    for attempt in range(1, fit_cfg.max_retries + 1):
        try:
            result = run_fit_once(
                image_data=image_data,
                adapter=adapter,
                specs=specs,
                sample_env=sample_env,
                fit_cfg=fit_cfg,
                rng=rng,
            )

            if result is None:
                raise RuntimeError("lmfit returned None")
            if not bool(getattr(result, "success", False)):
                raise RuntimeError(f"optimizer did not report success: {getattr(result, 'message', 'no message')}")
            if bool(getattr(result, "aborted", False)):
                raise RuntimeError("fit aborted")

            ok_unc, unc_msg = fit_result_has_valid_uncertainties(result, result.params)
            if not ok_unc:
                raise RuntimeError(unc_msg)

            fitted_params = maybe_relabel_two_sources(lmfit_params_to_dict(result.params))
            stderr_map = {name: float(result.params[name].stderr) for name in fitted_params.keys()}
            best_payload = {
                "sample_id": sample_id,
                "status": "success",
                "attempt": attempt,
                "message": str(getattr(result, "message", "")),
                "success": bool(getattr(result, "success", False)),
                "aborted": bool(getattr(result, "aborted", False)),
                "chisqr": float(getattr(result, "chisqr", np.nan)),
                "redchi": float(getattr(result, "redchi", np.nan)),
                "aic": float(getattr(result, "aic", np.nan)),
                "bic": float(getattr(result, "bic", np.nan)),
                "nfev": int(getattr(result, "nfev", -1)),
                "nvarys": int(getattr(result, "nvarys", -1)),
                "covar_ok": True,
                "errorbars_ok": True,
                "params": fitted_params,
                "stderr": stderr_map,
            }
            break
        except Exception as exc:  # noqa: PERF203
            last_exception = exc
            best_payload = {
                "sample_id": sample_id,
                "status": "retry_failed",
                "attempt": attempt,
                "message": str(exc),
                "success": False,
                "aborted": False,
                "chisqr": np.nan,
                "redchi": np.nan,
                "aic": np.nan,
                "bic": np.nan,
                "nfev": -1,
                "nvarys": -1,
                "covar_ok": False,
                "errorbars_ok": False,
                "params": {},
                "stderr": {},
                "traceback": traceback.format_exc(limit=2),
            }
    
    print(f"Sample {best_payload['sample_id']}, attempt {best_payload['attempt']}: fit completed with status={best_payload['status']}")
    return last_exception, best_payload



def fit_all(
    synthetic_manifest: pd.DataFrame,
    adapter: ModelAdapter,
    specs: Sequence[ParameterSpec],
    fit_cfg: FitSettings,
    pipeline_cfg: PipelineConfig,
    max_workers = 12,
) -> Tuple[pd.DataFrame, RunSummary]:
    
    base_dir = Path(pipeline_cfg.output_dir)
    dirs = ensure_dirs(base_dir)
    rng = np.random.default_rng(pipeline_cfg.random_seed + 1000)

    fit_rows: List[Dict[str, Any]] = []
    failed_ids: List[str] = []

    rec_list = [rec for _, rec in synthetic_manifest.iterrows()]
    workers = min(max_workers, os.cpu_count() or max_workers)

    # multiprocessing version 1
    '''
    with ProcessPoolExecutor(max_workers=workers) as executor:
        all_results = list(executor.map(run_fit_retrying, 
                                        rec_list, 
                                        [adapter] * len(rec_list), 
                                        [specs] * len(rec_list), 
                                        [fit_cfg] * len(rec_list), 
                                        [synthetic_manifest] * len(rec_list), 
                                        [rng] * len(rec_list)
                                        ))
    '''

    # multiprocessing version 2
    ray.init(num_cpus=max_workers, ignore_reinit_error=True)
    futures = [run_fit_retrying.remote(rec, adapter, specs, fit_cfg, synthetic_manifest, rng) for rec in rec_list]
    all_results = ray.get(futures)

    sample_id_list = [rec["sample_id"] for rec in rec_list]
    last_exception_list = [le for le, _ in all_results]
    best_payload_list = [bp for _, bp in all_results]
    result = {}
    for i in range(len(sample_id_list)):
        result[sample_id_list[i]] = last_exception_list[i], best_payload_list[i]

    for _, rec in synthetic_manifest.iterrows():
        sample_id = rec["sample_id"]
        # "for" cycle with attempts
        '''
        last_exception, best_payload = run_fit_retrying(rec, 
                                                        adapter, 
                                                        specs, 
                                                        fit_cfg, 
                                                        synthetic_manifest, 
                                                        rng
                                                        )
        '''
        last_exception, best_payload = result[sample_id]

        assert best_payload is not None

        if best_payload["status"] != "success":
            failed_ids.append(sample_id)
            best_payload["status"] = "failed_after_retries"
            if last_exception is not None:
                best_payload["message"] = str(last_exception)

        if pipeline_cfg.save_per_sample_fit_json:
            fit_json_path = dirs["fit_json"] / f"{sample_id}.json"
            fit_json_path.write_text(json.dumps(to_jsonable(best_payload), indent=2), encoding="utf-8")
        else:
            fit_json_path = None

        flat_row: Dict[str, Any] = {
            "sample_id": sample_id,
            "fit_status": best_payload["status"],
            "fit_attempt": best_payload["attempt"],
            "fit_message": best_payload["message"],
            "fit_json_path": str(fit_json_path) if fit_json_path else "",
            "chisqr": best_payload["chisqr"],
            "redchi": best_payload["redchi"],
            "aic": best_payload["aic"],
            "bic": best_payload["bic"],
            "nfev": best_payload["nfev"],
            "nvarys": best_payload["nvarys"],
            "covar_ok": best_payload["covar_ok"],
            "errorbars_ok": best_payload["errorbars_ok"],
        }
        for spec in specs:
            flat_row[f"fit_{spec.name}"] = best_payload["params"].get(spec.name, np.nan)
            flat_row[f"stderr_{spec.name}"] = best_payload["stderr"].get(spec.name, np.nan)
        fit_rows.append(flat_row)

    fit_df = pd.DataFrame(fit_rows)
    fit_df.to_csv(base_dir / "fit_manifest.csv", index=False)

    summary = RunSummary(
        n_total=int(len(fit_df)),
        n_success=int((fit_df["fit_status"] == "success").sum()),
        n_failed=int((fit_df["fit_status"] != "success").sum()),
        failed_sample_ids=failed_ids,
    )
    (base_dir / "fit_summary.json").write_text(json.dumps(to_jsonable(asdict(summary)), indent=2), encoding="utf-8")
    return fit_df, summary


# =============================================================================
# Блок 3. Сравнение truth vs fit и статистика
# =============================================================================


def robust_corr(a: np.ndarray, b: np.ndarray) -> Tuple[float, float]:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2:
        return np.nan, np.nan
    try:
        return float(pearsonr(a[mask], b[mask]).statistic), float(spearmanr(a[mask], b[mask]).statistic)
    except Exception:
        return np.nan, np.nan



def compute_parameter_statistics(df: pd.DataFrame, param_names: Sequence[str]) -> pd.DataFrame:
    rows = []
    for name in param_names:
        t = df[f"true_{name}"].to_numpy(dtype=float)
        f = df[f"fit_{name}"].to_numpy(dtype=float)
        s = df[f"stderr_{name}"].to_numpy(dtype=float)
        d = f - t
        mask = np.isfinite(t) & np.isfinite(f)
        pull_mask = mask & np.isfinite(s) & (s > 0)

        pearson_val, spearman_val = robust_corr(t, f)

        row = {
            "parameter": name,
            "n_valid": int(mask.sum()),
            "bias_mean": float(np.nanmean(d[mask])) if mask.any() else np.nan,
            "bias_median": float(np.nanmedian(d[mask])) if mask.any() else np.nan,
            "std_error": float(np.nanstd(d[mask], ddof=1)) if mask.sum() > 1 else np.nan,
            "mae": float(np.nanmean(np.abs(d[mask]))) if mask.any() else np.nan,
            "median_abs_error": float(np.nanmedian(np.abs(d[mask]))) if mask.any() else np.nan,
            "rmse": float(np.sqrt(np.nanmean(d[mask] ** 2))) if mask.any() else np.nan,
            "pearson_r": pearson_val,
            "spearman_r": spearman_val,
            "mean_true": float(np.nanmean(t[mask])) if mask.any() else np.nan,
            "mean_fit": float(np.nanmean(f[mask])) if mask.any() else np.nan,
            "std_true": float(np.nanstd(t[mask], ddof=1)) if mask.sum() > 1 else np.nan,
            "std_fit": float(np.nanstd(f[mask], ddof=1)) if mask.sum() > 1 else np.nan,
            "coverage_1sigma": float(np.nanmean(np.abs(d[pull_mask]) <= s[pull_mask])) if pull_mask.any() else np.nan,
            "coverage_2sigma": float(np.nanmean(np.abs(d[pull_mask]) <= 2.0 * s[pull_mask])) if pull_mask.any() else np.nan,
            "pull_mean": float(np.nanmean((d[pull_mask] / s[pull_mask]))) if pull_mask.any() else np.nan,
            "pull_std": float(np.nanstd((d[pull_mask] / s[pull_mask]), ddof=1)) if pull_mask.sum() > 1 else np.nan,
        }
        rows.append(row)
    return pd.DataFrame(rows)



def save_histograms_overlay(
    df: pd.DataFrame,
    param_names: Sequence[str],
    out_path: Path,
    prefix_true: str,
    prefix_fit: str,
) -> None:
    n = len(param_names)
    ncols = 3
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.8 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, param_names):
        ax.hist(df[f"{prefix_true}{name}"].dropna().to_numpy(), bins=20, alpha=0.55, label="true")
        ax.hist(df[f"{prefix_fit}{name}"].dropna().to_numpy(), bins=20, alpha=0.55, label="fit")
        ax.set_title(name)
        ax.legend()
        ax.grid(alpha=0.25)
    for ax in axes[len(param_names):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)



def save_error_histograms(df: pd.DataFrame, param_names: Sequence[str], out_path: Path) -> None:
    n = len(param_names)
    ncols = 3
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.8 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, param_names):
        delta = df[f"fit_{name}"] - df[f"true_{name}"]
        ax.hist(delta.dropna().to_numpy(), bins=20, alpha=0.75)
        ax.axvline(0.0, linestyle="--", linewidth=1.0)
        ax.set_title(f"Δ {name}")
        ax.grid(alpha=0.25)
    for ax in axes[len(param_names):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)



def save_pull_histograms(df: pd.DataFrame, param_names: Sequence[str], out_path: Path) -> None:
    n = len(param_names)
    ncols = 3
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.8 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, param_names):
        pull = (df[f"fit_{name}"] - df[f"true_{name}"]) / df[f"stderr_{name}"]
        pull = pull.replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
        ax.hist(pull, bins=20, alpha=0.75)
        ax.axvline(0.0, linestyle="--", linewidth=1.0)
        ax.set_title(f"pull {name}")
        ax.grid(alpha=0.25)
    for ax in axes[len(param_names):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)



def save_true_vs_fit_scatter(df: pd.DataFrame, param_names: Sequence[str], out_path: Path) -> None:
    n = len(param_names)
    ncols = 3
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.2 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, param_names):
        x = df[f"true_{name}"].to_numpy(dtype=float)
        y = df[f"fit_{name}"].to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.any():
            ax.scatter(x[mask], y[mask], s=14, alpha=0.65)
            lo = float(np.nanmin(np.r_[x[mask], y[mask]]))
            hi = float(np.nanmax(np.r_[x[mask], y[mask]]))
            ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0)
        ax.set_xlabel(f"true {name}")
        ax.set_ylabel(f"fit {name}")
        ax.set_title(name)
        ax.grid(alpha=0.25)
    for ax in axes[len(param_names):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)



def save_goodness_histograms(df: pd.DataFrame, out_path: Path) -> None:
    metrics = ["chisqr", "redchi", "aic", "bic"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes = axes.ravel()
    for ax, metric in zip(axes, metrics):
        vals = df[metric].replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
        ax.hist(vals, bins=20, alpha=0.75)
        ax.set_title(metric)
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)



def save_error_correlation_heatmap(df: pd.DataFrame, param_names: Sequence[str], out_path: Path) -> None:
    err_cols = [f"delta_{name}" for name in param_names]
    err_df = pd.DataFrame({col: df[col] for col in err_cols})
    corr = err_df.corr(numeric_only=True)
    fig, ax = plt.subplots(figsize=(0.7 * len(param_names) + 4, 0.7 * len(param_names) + 4))
    im = ax.imshow(corr.to_numpy(dtype=float), aspect="auto")
    ax.set_xticks(np.arange(len(param_names)))
    ax.set_yticks(np.arange(len(param_names)))
    ax.set_xticklabels(param_names, rotation=90)
    ax.set_yticklabels(param_names)
    fig.colorbar(im, ax=ax, shrink=0.85)
    ax.set_title("Correlation matrix of parameter errors")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)



def save_example_residual_plots(
    merged: pd.DataFrame,
    adapter: ModelAdapter,
    param_names: Sequence[str],
    out_dir: Path,
    n_examples: int,
) -> None:
    successful = merged.loc[merged["fit_status"] == "success"].head(n_examples)
    for _, row in successful.iterrows():
        sample_id = row["sample_id"]
        image = np.load(row["image_path"]).astype(np.float64)
        env = {col.replace("env_", ""): row[col] for col in merged.columns if col.startswith("env_")}
        fit_env = dict(env)
        fit_env["noise_sigma"] = 0.0
        fit_params = {name: float(row[f"fit_{name}"]) for name in param_names}
        model = adapter(fit_params, fit_env, rng=None)
        resid = image - model

        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
        for ax, arr, title in zip(axes, [image, model, resid], ["data", "best-fit model", "residual"]):
            im = ax.imshow(arr, origin="lower")
            ax.set_title(f"{sample_id}: {title}")
            fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        fig.savefig(out_dir / f"{sample_id}_triptych.png", dpi=160)
        plt.close(fig)


def resolve_data_path(path_like: Any, base_dir: Path) -> Path:
    path = Path(str(path_like))
    if path.is_absolute():
        return path
    if path.exists():
        return path
    candidate = base_dir / path
    if candidate.exists():
        return candidate
    return path


def save_fit_json_image_plots(
    synthetic_manifest: pd.DataFrame,
    fit_manifest: pd.DataFrame,
    adapter: ModelAdapter,
    param_names: Sequence[str],
    pipeline_cfg: PipelineConfig,
) -> int:
    """Save data/model/residual plots using fit JSON params, without truth columns."""
    base_dir = Path(pipeline_cfg.output_dir)
    dirs = ensure_dirs(base_dir)
    manifest_by_sample = synthetic_manifest.set_index("sample_id", drop=False)
    saved_count = 0

    for _, fit_row in fit_manifest.iterrows():
        sample_id = fit_row["sample_id"]
        if fit_row.get("fit_status", "") != "success":
            continue
        if sample_id not in manifest_by_sample.index:
            print(f"[makeimg] skipping {sample_id}: not found in synthetic_manifest.csv")
            continue

        data_row = manifest_by_sample.loc[sample_id]
        image_path = resolve_data_path(data_row["image_path"], base_dir)
        fit_json_path = resolve_data_path(fit_row.get("fit_json_path", ""), base_dir)
        if not image_path.exists():
            print(f"[makeimg] skipping {sample_id}: image not found: {image_path}")
            continue
        if not fit_json_path.exists():
            print(f"[makeimg] skipping {sample_id}: fit JSON not found: {fit_json_path}")
            continue

        payload = json.loads(fit_json_path.read_text(encoding="utf-8"))
        fit_params = {name: float(payload["params"][name]) for name in param_names}
        image = np.load(image_path).astype(np.float64)
        env = {col.replace("env_", ""): data_row[col] for col in synthetic_manifest.columns if col.startswith("env_")}
        fit_env = dict(env)
        fit_env["noise_sigma"] = 0.0
        model = adapter(fit_params, fit_env, rng=None)
        resid = image - model

        fig, axes = plt.subplots(3, 1, figsize=(6.4, 12.0))
        for ax, arr, title in zip(axes, [image, model, resid], ["input image", "best-fit model", "residual"]):
            im = ax.imshow(arr, origin="lower")
            ax.set_title(f"{sample_id}: {title}")
            fig.colorbar(im, ax=ax, shrink=0.82)
        fig.tight_layout()
        fig.savefig(dirs["examples"] / f"{sample_id}_fit_images.png", dpi=160)
        plt.close(fig)
        saved_count += 1

    print(f"[makeimg] saved {saved_count} image plot(s) to {dirs['examples']}")
    return saved_count



def analyze_results(
    synthetic_manifest: pd.DataFrame,
    fit_manifest: pd.DataFrame,
    specs: Sequence[ParameterSpec],
    adapter: ModelAdapter,
    pipeline_cfg: PipelineConfig,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    base_dir = Path(pipeline_cfg.output_dir)
    dirs = ensure_dirs(base_dir)
    param_names = [spec.name for spec in specs]

    merged = synthetic_manifest.merge(fit_manifest, on="sample_id", how="left")
    merged_success = merged.loc[merged["fit_status"] == "success"].copy()

    for name in param_names:
        merged_success[f"delta_{name}"] = merged_success[f"fit_{name}"] - merged_success[f"true_{name}"]
        merged_success[f"pull_{name}"] = merged_success[f"delta_{name}"] / merged_success[f"stderr_{name}"]

    parameter_stats = compute_parameter_statistics(merged_success, param_names)
    parameter_stats.to_csv(dirs["tables"] / "parameter_statistics.csv", index=False)
    merged.to_csv(dirs["tables"] / "merged_truth_fit.csv", index=False)

    global_stats = {
        "n_total": int(len(merged)),
        "n_success": int((merged["fit_status"] == "success").sum()),
        "n_failed": int((merged["fit_status"] != "success").sum()),
        "success_fraction": float((merged["fit_status"] == "success").mean()) if len(merged) else np.nan,
        "median_redchi": float(np.nanmedian(merged_success["redchi"])) if len(merged_success) else np.nan,
        "mean_redchi": float(np.nanmean(merged_success["redchi"])) if len(merged_success) else np.nan,
        "median_chisqr": float(np.nanmedian(merged_success["chisqr"])) if len(merged_success) else np.nan,
        "median_aic": float(np.nanmedian(merged_success["aic"])) if len(merged_success) else np.nan,
        "median_bic": float(np.nanmedian(merged_success["bic"])) if len(merged_success) else np.nan,
        "median_attempts": float(np.nanmedian(merged.loc[merged["fit_status"] == "success", "fit_attempt"])) if len(merged_success) else np.nan,
    }
    (dirs["tables"] / "global_statistics.json").write_text(json.dumps(to_jsonable(global_stats), indent=2), encoding="utf-8")

    if len(merged_success):
        save_histograms_overlay(merged_success, param_names, dirs["plots"] / "hist_true_vs_fit.png", "true_", "fit_")
        save_error_histograms(merged_success, param_names, dirs["plots"] / "hist_parameter_errors.png")
        save_pull_histograms(merged_success, param_names, dirs["plots"] / "hist_parameter_pulls.png")
        save_true_vs_fit_scatter(merged_success, param_names, dirs["plots"] / "scatter_true_vs_fit.png")
        save_goodness_histograms(merged_success, dirs["plots"] / "goodness_of_fit_histograms.png")
        save_error_correlation_heatmap(merged_success, param_names, dirs["plots"] / "error_correlation_heatmap.png")
        if pipeline_cfg.save_example_images:
            save_example_residual_plots(
                merged=merged,
                adapter=adapter,
                param_names=param_names,
                out_dir=dirs["examples"],
                n_examples=pipeline_cfg.example_plots_to_save,
            )

    failed_df = merged.loc[merged["fit_status"] != "success", ["sample_id", "fit_status", "fit_attempt", "fit_message"]].copy()
    failed_df.to_csv(dirs["tables"] / "failed_samples.csv", index=False)

    report_lines = [
        "# Fit quality report",
        f"Total synthetic images: {global_stats['n_total']}",
        f"Successful fits: {global_stats['n_success']}",
        f"Failed fits: {global_stats['n_failed']}",
        f"Success fraction: {global_stats['success_fraction']:.4f}" if np.isfinite(global_stats['success_fraction']) else "Success fraction: nan",
        f"Median reduced chi-square: {global_stats['median_redchi']:.6g}" if np.isfinite(global_stats['median_redchi']) else "Median reduced chi-square: nan",
        "",
        "Saved files:",
        "- synthetic_manifest.csv",
        "- fit_manifest.csv",
        "- analysis/tables/parameter_statistics.csv",
        "- analysis/tables/global_statistics.json",
        "- analysis/tables/failed_samples.csv",
        "- analysis/plots/*.png",
    ]
    (base_dir / "analysis_report.md").write_text("\n".join(report_lines), encoding="utf-8")
    return parameter_stats, global_stats


# =============================================================================
# Сценарий запуска
# =============================================================================


def build_default_environment() -> EnvDict:
    return {
        "nx": 128,
        "ny": 128,
        "pixel_scale": 0.10,   # arcsec/pixel
        "noise_sigma": 0.1,   # ненулевой шум
        "psf_sigma": 0.20,     # arcsec, ненулевая ширина PSF
    }



def run_pipeline(
    adapter: ModelAdapter,
    specs: Sequence[ParameterSpec],
    env: EnvDict,
    pipeline_cfg: PipelineConfig,
    fit_cfg: FitSettings,
    mode: str = "all",
    makeimg: bool = False,
) -> None:
    synthetic_manifest_path = Path(pipeline_cfg.output_dir) / "synthetic_manifest.csv"
    fit_manifest_path = Path(pipeline_cfg.output_dir) / "fit_manifest.csv"

    synthetic_manifest = None
    fit_manifest = None

    if mode in {"generate", "all"}:
        synthetic_manifest = generate_synthetic_dataset(adapter, specs, env, pipeline_cfg)
    elif synthetic_manifest_path.exists():
        synthetic_manifest = pd.read_csv(synthetic_manifest_path)
    else:
        raise FileNotFoundError("synthetic_manifest.csv not found. Run with --mode generate or --mode all first.")

    if mode in {"fit", "all"}:
        fit_manifest, summary = fit_all(synthetic_manifest, adapter, specs, fit_cfg, pipeline_cfg)
        print(f"[fit] success={summary.n_success}, failed={summary.n_failed}")
    elif fit_manifest_path.exists():
        fit_manifest = pd.read_csv(fit_manifest_path)
    else:
        raise FileNotFoundError("fit_manifest.csv not found. Run with --mode fit or --mode all first.")

    if mode in {"analyze", "all"}:
        parameter_stats, global_stats = analyze_results(synthetic_manifest, fit_manifest, specs, adapter, pipeline_cfg)
        print("[analyze] global stats:")
        print(json.dumps(to_jsonable(global_stats), indent=2))
        print("[analyze] first rows of parameter statistics:")
        print(parameter_stats.head())

    if makeimg:
        save_fit_json_image_plots(
            synthetic_manifest=synthetic_manifest,
            fit_manifest=fit_manifest,
            adapter=adapter,
            param_names=[spec.name for spec in specs],
            pipeline_cfg=pipeline_cfg,
        )



def parse_args() -> argparse.Namespace:
    env = build_default_environment()
    parser = argparse.ArgumentParser(description="Pipeline for fit-quality validation on synthetic 2D images.")
    parser.add_argument("--mode", choices=["generate", "fit", "analyze", "all", "makeimg", "real"], default="all")
    parser.add_argument("--output-dir", default="SL_result")
    parser.add_argument("--real-input-dir", help="Prepared dataset from photometry_cli.py for --mode real")
    parser.add_argument("--real-output-dir", default="../result_real_pointfit")
    parser.add_argument("--real-starts", type=int, default=6)
    parser.add_argument("--real-rounds", type=int, default=2)
    parser.add_argument("--n-samples", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-retries", type=int, default=7)
    parser.add_argument("--methods", nargs="+", default=["least_squares", "nelder", "leastsq"], help="Chain of lmfit methods.")
    parser.add_argument("--max-nfev", type=int, default=15000)
    parser.add_argument("--nx", type=int, default=env["nx"])
    parser.add_argument("--ny", type=int, default=env["ny"])
    parser.add_argument("--pixel-scale", type=float, default=env["pixel_scale"])
    parser.add_argument("--noise-sigma", type=float, default=env["noise_sigma"])
    parser.add_argument("--psf-sigma", type=float, default=env["psf_sigma"])
    parser.add_argument("--example-plots", type=int, default=10)
    parser.add_argument(
        "--makeimg",
        action="store_true",
        help="Save input/model/residual PNGs from fit_results/json into analysis/examples.",
    )
    return parser.parse_args()



def main() -> None:
    args = parse_args()
    mode = args.mode
    if args.makeimg and "--mode" not in sys.argv:
        mode = "makeimg"

    env = build_default_environment()
    env.update(
        {
            "nx": args.nx,
            "ny": args.ny,
            "pixel_scale": args.pixel_scale,
            "noise_sigma": args.noise_sigma,
            "psf_sigma": args.psf_sigma,
        }
    )

    pipeline_cfg = PipelineConfig(
        output_dir=args.output_dir,
        n_samples=args.n_samples,
        random_seed=args.seed,
        example_plots_to_save=args.example_plots,
    )
    fit_cfg = FitSettings(
        methods_chain=args.methods,
        max_retries=args.max_retries,
        max_nfev=args.max_nfev,
    )

    specs = build_real_parameter_specs(env)
    adapter = ModelAdapter(
        external_func=real_external_galaxy_plus_string,
        model_param_names=[spec.name for spec in specs],
        env_param_names=["nx", "ny", "pixel_scale", "noise_sigma", "psf_sigma"],
    )

    if mode == "real":
        if not args.real_input_dir:
            raise SystemExit("--mode real требует --real-input-dir")
        subprocess.run(
            [sys.executable, str(Path(__file__).with_name("real_fit.py")),
             "--input-dir", args.real_input_dir,
             "--output-dir", args.real_output_dir,
             "--starts", str(args.real_starts),
             "--rounds", str(args.real_rounds),
             "--seed", str(args.seed)],
            check=True,
        )
        return

    run_pipeline(
        adapter=adapter,
        specs=specs,
        env=env,
        pipeline_cfg=pipeline_cfg,
        fit_cfg=fit_cfg,
        mode=mode,
        makeimg=args.makeimg or mode == "makeimg",
    )


if __name__ == "__main__":
    '''
    env = build_default_environment()
    img_test = real_external_galaxy_plus_string(nx = env["nx"], 
                                                ny = env["ny"], 
                                                x0_pix = 0.0,
                                                y0_pix = 0.0,
                                                pos_angle_rad = 0.0,
                                                pos_angle_string_rad = 0.0,
                                                pixel_scale = env["pixel_scale"],
                                                noise_sigma = env["noise_sigma"],
                                                psf_sigma = env["psf_sigma"],
                                                )
    string_params = {
        'tension': 1e-6,
        'inclination': 0.0,
        'pos_angle_string_rad': 0.0,
        'distane_center_string': 0.0,
        'RsRg': 0.5
    }
    sl_utils.plot_simple_image(img_test, env["pixel_scale"], False, string_params)
    '''
    main()
