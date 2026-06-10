from __future__ import annotations

import argparse
import json
import math
import os
import traceback
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib-telescope-simulation"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fit_quality_pipeline as fqp


IE_MIN = 1.0
IE_MAX = 10.0**1.6
NEAR_TRUE_REL_TOL = 0.20
NEAR_ZERO_RANGE_FRACTION = 0.05

TELESCOPE_CONFIGS = {
    "DOT": "params_DOT.json",
    "KECK": "params_KECK.json",
}

MAIN_STRING_PARAMETERS = [
    "ThetaE",
    "dThetaEdXi",
    "pos_angle_string_rad",
    "distane_center_string",
]

PREFERRED_PSF_CASES = [
    "site_median_seeing",
    "combined_seeing_0p7_plus_instrumental_center",
    "seeing_0p7_arcsec",
    "median_observed_V_2016_2021",
]


def uniform_sampler(lo: float, hi: float):
    return lambda rng, _env: rng.uniform(lo, hi)


def sky_noisy_real_external_galaxy_plus_string(
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
) -> np.ndarray:
    image = fqp.real_external_galaxy_plus_string(
        ny,
        nx,
        Ie=Ie,
        re_arcsec=re_arcsec,
        n=n,
        x0_pix=x0_pix,
        y0_pix=y0_pix,
        q=q,
        pos_angle_rad=pos_angle_rad,
        ThetaE=ThetaE,
        dThetaEdXi=dThetaEdXi,
        pos_angle_string_rad=pos_angle_string_rad,
        distane_center_string=distane_center_string,
        pixel_scale=pixel_scale,
        noise_sigma=0.0,
        psf_sigma=psf_sigma,
        rng=None,
    )
    if noise_sigma > 0:
        if rng is None:
            rng = np.random.default_rng()
        image = image + rng.normal(loc=0.0, scale=float(noise_sigma), size=image.shape)
    return np.asarray(image, dtype=np.float64)


def build_telescope_parameter_specs(env: fqp.EnvDict) -> List[fqp.ParameterSpec]:
    specs = fqp.build_real_parameter_specs(env)
    return [
        fqp.ParameterSpec(
            "Ie",
            IE_MIN,
            IE_MAX,
            uniform_sampler(IE_MIN, IE_MAX),
        )
        if spec.name == "Ie"
        else spec
        for spec in specs
    ]


def build_adapter(specs: Sequence[fqp.ParameterSpec]) -> fqp.ModelAdapter:
    return fqp.ModelAdapter(
        external_func=sky_noisy_real_external_galaxy_plus_string,
        model_param_names=[spec.name for spec in specs],
        env_param_names=["nx", "ny", "pixel_scale", "noise_sigma", "psf_sigma"],
    )


def load_telescope_config(repo_root: Path, telescope: str) -> Dict[str, Any]:
    path = repo_root / TELESCOPE_CONFIGS[telescope]
    if not path.exists():
        raise FileNotFoundError(f"Telescope parameter file not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def valid_sky_filters(config: Dict[str, Any]) -> Dict[str, float]:
    filters = config.get("sky_background", {}).get("filters", {})
    out: Dict[str, float] = {}
    for name, payload in filters.items():
        rate = payload.get("sky_electrons_per_second_per_pixel")
        if isinstance(rate, (int, float)) and np.isfinite(rate):
            out[name] = float(rate)
    return out


def exposure_grid(t_exp_min: float, t_exp_max: float, n_exposures: int) -> np.ndarray:
    if t_exp_min <= 0 or t_exp_max <= 0:
        raise ValueError("Exposure times must be positive.")
    if n_exposures <= 0:
        raise ValueError("--n-exposures must be positive.")
    if n_exposures == 1:
        return np.asarray([float(t_exp_min)], dtype=float)
    return 10.0 ** np.linspace(math.log10(t_exp_min), math.log10(t_exp_max), n_exposures)


def exposure_dir_name(t_exp: float) -> str:
    text = f"{float(t_exp):.4g}"
    if "e" in text.lower():
        text = f"{float(t_exp):.6f}".rstrip("0").rstrip(".")
    text = text.replace(".", "p").replace("+", "").replace("-", "m")
    return f"t_exp_{text}s"


def sky_sigma_image_units(sky_rate_e_per_s_pix: float, t_exp: float) -> float:
    return float(math.sqrt(sky_rate_e_per_s_pix * t_exp) / t_exp)


def choose_psf_sigma(config: Dict[str, Any], psf_case: Optional[str], psf_sigma_override: Optional[float]) -> float:
    if psf_sigma_override is not None:
        return float(psf_sigma_override)

    cases = config.get("psf", {}).get("cases", {})
    if not cases:
        return float(fqp.build_default_environment()["psf_sigma"])

    if psf_case is not None:
        if psf_case not in cases:
            available = ", ".join(cases.keys())
            raise ValueError(f"PSF case '{psf_case}' is not in config. Available cases: {available}")
        return float(cases[psf_case]["sigma_arcsec"])

    for case_name in PREFERRED_PSF_CASES:
        if case_name in cases and cases[case_name].get("sigma_arcsec") is not None:
            return float(cases[case_name]["sigma_arcsec"])

    for payload in cases.values():
        if payload.get("sigma_arcsec") is not None:
            return float(payload["sigma_arcsec"])

    return float(fqp.build_default_environment()["psf_sigma"])


def build_environment(
    config: Dict[str, Any],
    sky_rate: float,
    t_exp: float,
    nx: int,
    ny: int,
    psf_case: Optional[str],
    psf_sigma_override: Optional[float],
) -> fqp.EnvDict:
    env = fqp.build_default_environment()
    env.update(
        {
            "nx": int(nx),
            "ny": int(ny),
            "pixel_scale": float(config["instrument"]["pixel_scale_arcsec_per_pixel"]),
            "noise_sigma": sky_sigma_image_units(sky_rate, t_exp),
            "psf_sigma": choose_psf_sigma(config, psf_case, psf_sigma_override),
            "t_exp": float(t_exp),
            "sky_electrons_per_second_per_pixel": float(sky_rate),
        }
    )
    return env


def write_run_metadata(
    run_dir: Path,
    telescope: str,
    filter_name: str,
    t_exp: float,
    sky_rate: float,
    env: fqp.EnvDict,
) -> None:
    payload = {
        "telescope": telescope,
        "filter": filter_name,
        "t_exp_seconds": float(t_exp),
        "sky_electrons_per_second_per_pixel": float(sky_rate),
        "sky_variance_electrons_per_pixel": float(sky_rate * t_exp),
        "sky_sigma_electrons_per_pixel": float(math.sqrt(sky_rate * t_exp)),
        "noise_sigma_image_units": float(env["noise_sigma"]),
        "Ie_min": IE_MIN,
        "Ie_max": IE_MAX,
        "env": env,
    }
    (run_dir / "run_metadata.json").write_text(
        json.dumps(fqp.to_jsonable(payload), indent=2),
        encoding="utf-8",
    )


def run_fit_retrying_serial(
    rec: pd.Series,
    adapter: fqp.ModelAdapter,
    specs: Sequence[fqp.ParameterSpec],
    fit_cfg: fqp.FitSettings,
    synthetic_manifest: pd.DataFrame,
    rng: np.random.Generator,
) -> Tuple[Optional[Exception], Dict[str, Any]]:
    sample_id = rec["sample_id"]
    image_data = np.load(rec["image_path"]).astype(np.float64)
    sample_env = {
        col.replace("env_", ""): rec[col]
        for col in synthetic_manifest.columns
        if col.startswith("env_")
    }

    best_payload: Optional[Dict[str, Any]] = None
    last_exception: Optional[Exception] = None

    for attempt in range(1, fit_cfg.max_retries + 1):
        try:
            result = fqp.run_fit_once(
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

            ok_unc, unc_msg = fqp.fit_result_has_valid_uncertainties(result, result.params)
            if not ok_unc:
                raise RuntimeError(unc_msg)

            fitted_params = fqp.maybe_relabel_two_sources(fqp.lmfit_params_to_dict(result.params))
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
        except Exception as exc:
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

    assert best_payload is not None
    print(f"Sample {best_payload['sample_id']}, attempt {best_payload['attempt']}: fit completed with status={best_payload['status']}")
    return last_exception, best_payload


def fit_all_serial(
    synthetic_manifest: pd.DataFrame,
    adapter: fqp.ModelAdapter,
    specs: Sequence[fqp.ParameterSpec],
    fit_cfg: fqp.FitSettings,
    pipeline_cfg: fqp.PipelineConfig,
) -> Tuple[pd.DataFrame, fqp.RunSummary]:
    base_dir = Path(pipeline_cfg.output_dir)
    dirs = fqp.ensure_dirs(base_dir)
    rng = np.random.default_rng(pipeline_cfg.random_seed + 1000)
    fit_rows: List[Dict[str, Any]] = []
    failed_ids: List[str] = []

    for _, rec in synthetic_manifest.iterrows():
        sample_id = rec["sample_id"]
        last_exception, best_payload = run_fit_retrying_serial(
            rec=rec,
            adapter=adapter,
            specs=specs,
            fit_cfg=fit_cfg,
            synthetic_manifest=synthetic_manifest,
            rng=rng,
        )

        if best_payload["status"] != "success":
            failed_ids.append(sample_id)
            best_payload["status"] = "failed_after_retries"
            if last_exception is not None:
                best_payload["message"] = str(last_exception)

        if pipeline_cfg.save_per_sample_fit_json:
            fit_json_path = dirs["fit_json"] / f"{sample_id}.json"
            fit_json_path.write_text(json.dumps(fqp.to_jsonable(best_payload), indent=2), encoding="utf-8")
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

    summary = fqp.RunSummary(
        n_total=int(len(fit_df)),
        n_success=int((fit_df["fit_status"] == "success").sum()),
        n_failed=int((fit_df["fit_status"] != "success").sum()),
        failed_sample_ids=failed_ids,
    )
    (base_dir / "fit_summary.json").write_text(
        json.dumps(fqp.to_jsonable(asdict(summary)), indent=2),
        encoding="utf-8",
    )
    return fit_df, summary


def run_pipeline_point(
    run_dir: Path,
    env: fqp.EnvDict,
    n_samples: int,
    seed: int,
    example_plots: int,
    max_retries: int,
    methods: Sequence[str],
    max_nfev: int,
    max_workers: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[fqp.ParameterSpec]]:
    specs = build_telescope_parameter_specs(env)
    adapter = build_adapter(specs)
    pipeline_cfg = fqp.PipelineConfig(
        output_dir=str(run_dir),
        n_samples=int(n_samples),
        random_seed=int(seed),
        example_plots_to_save=int(example_plots),
        save_per_sample_fit_json=True,
        save_example_images=True,
    )
    fit_cfg = fqp.FitSettings(
        methods_chain=list(methods),
        max_retries=int(max_retries),
        max_nfev=int(max_nfev),
    )

    synthetic_manifest = fqp.generate_synthetic_dataset(adapter, specs, env, pipeline_cfg)
    if max_workers <= 1:
        fit_manifest, summary = fit_all_serial(
            synthetic_manifest=synthetic_manifest,
            adapter=adapter,
            specs=specs,
            fit_cfg=fit_cfg,
            pipeline_cfg=pipeline_cfg,
        )
    else:
        fit_manifest, summary = fqp.fit_all(
            synthetic_manifest=synthetic_manifest,
            adapter=adapter,
            specs=specs,
            fit_cfg=fit_cfg,
            pipeline_cfg=pipeline_cfg,
            max_workers=int(max_workers),
        )
        fqp.ray.shutdown()
    print(f"[fit] {run_dir}: success={summary.n_success}, failed={summary.n_failed}")
    fqp.analyze_results(
        synthetic_manifest=synthetic_manifest,
        fit_manifest=fit_manifest,
        specs=specs,
        adapter=adapter,
        pipeline_cfg=pipeline_cfg,
    )
    return synthetic_manifest, fit_manifest, specs


def parameter_range(specs: Sequence[fqp.ParameterSpec], name: str) -> float:
    for spec in specs:
        if spec.name == name:
            return float(spec.max_value - spec.min_value)
    raise KeyError(name)


def wrapped_parameter_error(name: str, fit_value: np.ndarray, true_value: np.ndarray) -> np.ndarray:
    diff = np.abs(fit_value - true_value)
    if name == "pos_angle_string_rad":
        period = np.pi
        diff = np.abs((fit_value - true_value + period / 2.0) % period - period / 2.0)
    return diff


def recovery_frame(merged: pd.DataFrame, specs: Sequence[fqp.ParameterSpec]) -> pd.DataFrame:
    out = merged.copy()
    recoveries = []
    near_true_mask = out["fit_status"].eq("success").to_numpy().copy()

    for name in MAIN_STRING_PARAMETERS:
        true_value = out[f"true_{name}"].to_numpy(dtype=float)
        fit_value = out[f"fit_{name}"].to_numpy(dtype=float)
        abs_error = wrapped_parameter_error(name, fit_value, true_value)
        fallback = NEAR_ZERO_RANGE_FRACTION * parameter_range(specs, name)
        denominator = np.maximum(np.abs(true_value), fallback)
        rel_error = abs_error / denominator
        recovery = 100.0 * (1.0 - rel_error)
        out[f"rel_error_{name}"] = rel_error
        out[f"recovery_{name}"] = recovery
        recoveries.append(recovery)
        near_true_mask &= np.isfinite(rel_error) & (rel_error <= NEAR_TRUE_REL_TOL)

    out["main_string_recovery_percent"] = np.nanmean(np.vstack(recoveries), axis=0)
    out["near_true"] = near_true_mask
    return out


def mean_and_sem(values: Iterable[float]) -> Tuple[float, float]:
    arr = pd.Series(values, dtype="float64").replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    if len(arr) == 0:
        return np.nan, np.nan
    if len(arr) == 1:
        return float(arr[0]), np.nan
    return float(np.mean(arr)), float(np.std(arr, ddof=1) / math.sqrt(len(arr)))


def summarize_point(
    synthetic_manifest: pd.DataFrame,
    fit_manifest: pd.DataFrame,
    specs: Sequence[fqp.ParameterSpec],
    telescope: str,
    filter_name: str,
    t_exp: float,
    sky_rate: float,
    noise_sigma: float,
    run_dir: Path,
) -> Dict[str, Any]:
    merged = synthetic_manifest.merge(fit_manifest, on="sample_id", how="left")
    merged = recovery_frame(merged, specs)
    success = merged["fit_status"].eq("success")
    near_true = success & merged["near_true"]

    redchi_success_mean, redchi_success_sem = mean_and_sem(merged.loc[success, "redchi"])
    redchi_near_mean, redchi_near_sem = mean_and_sem(merged.loc[near_true, "redchi"])
    rec_success_mean, rec_success_sem = mean_and_sem(merged.loc[success, "main_string_recovery_percent"])
    rec_near_mean, rec_near_sem = mean_and_sem(merged.loc[near_true, "main_string_recovery_percent"])

    return {
        "telescope": telescope,
        "filter": filter_name,
        "t_exp_seconds": float(t_exp),
        "run_dir": str(run_dir),
        "sky_electrons_per_second_per_pixel": float(sky_rate),
        "noise_sigma_image_units": float(noise_sigma),
        "n_total": int(len(merged)),
        "n_success": int(success.sum()),
        "n_failed": int((~success).sum()),
        "success_percent": float(100.0 * success.mean()) if len(merged) else np.nan,
        "n_near_true": int(near_true.sum()),
        "near_true_percent": float(100.0 * near_true.mean()) if len(merged) else np.nan,
        "redchi_success_mean": redchi_success_mean,
        "redchi_success_sem": redchi_success_sem,
        "redchi_near_true_mean": redchi_near_mean,
        "redchi_near_true_sem": redchi_near_sem,
        "recovery_success_mean_percent": rec_success_mean,
        "recovery_success_sem_percent": rec_success_sem,
        "recovery_near_true_mean_percent": rec_near_mean,
        "recovery_near_true_sem_percent": rec_near_sem,
    }


def save_errorbar_plot(
    df: pd.DataFrame,
    y_col: str,
    yerr_col: Optional[str],
    out_path: Path,
    title: str,
    ylabel: str,
    compare_col: Optional[str] = None,
    compare_yerr_col: Optional[str] = None,
    compare_label: str = "near-true",
) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    x = df["t_exp_seconds"].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)
    yerr = df[yerr_col].to_numpy(dtype=float) if yerr_col else None
    ax.errorbar(x, y, yerr=yerr, marker="o", linewidth=1.6, capsize=3, label="successful")

    if compare_col is not None:
        y_compare = df[compare_col].to_numpy(dtype=float)
        yerr_compare = df[compare_yerr_col].to_numpy(dtype=float) if compare_yerr_col else None
        ax.errorbar(
            x,
            y_compare,
            yerr=yerr_compare,
            marker="s",
            linewidth=1.6,
            capsize=3,
            label=compare_label,
        )
        ax.legend()

    ax.set_xscale("log")
    ax.set_xlabel("Exposure time, s")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def save_filter_results(summary_df: pd.DataFrame, results_dir: Path, telescope: str, filter_name: str) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    summary_df = summary_df.sort_values("t_exp_seconds").reset_index(drop=True)
    summary_df.to_csv(results_dir / "simulation_summary.csv", index=False)
    (results_dir / "simulation_summary.json").write_text(
        json.dumps(fqp.to_jsonable(summary_df.to_dict(orient="records")), indent=2),
        encoding="utf-8",
    )

    prefix = f"{telescope} {filter_name}"
    save_errorbar_plot(
        summary_df,
        "success_percent",
        None,
        results_dir / "success_percent_vs_t_exp.png",
        f"{prefix}: successful fits",
        "Successful fits, %",
    )
    save_errorbar_plot(
        summary_df,
        "redchi_success_mean",
        "redchi_success_sem",
        results_dir / "redchi_success_vs_t_exp.png",
        f"{prefix}: reduced chi-square",
        "Mean reduced chi-square",
    )
    save_errorbar_plot(
        summary_df,
        "recovery_success_mean_percent",
        "recovery_success_sem_percent",
        results_dir / "recovery_success_vs_t_exp.png",
        f"{prefix}: main string recovery",
        "Mean recovery, %",
    )
    save_errorbar_plot(
        summary_df,
        "redchi_success_mean",
        "redchi_success_sem",
        results_dir / "redchi_success_and_near_true_vs_t_exp.png",
        f"{prefix}: reduced chi-square",
        "Mean reduced chi-square",
        compare_col="redchi_near_true_mean",
        compare_yerr_col="redchi_near_true_sem",
    )
    save_errorbar_plot(
        summary_df,
        "recovery_success_mean_percent",
        "recovery_success_sem_percent",
        results_dir / "recovery_success_and_near_true_vs_t_exp.png",
        f"{prefix}: main string recovery",
        "Mean recovery, %",
        compare_col="recovery_near_true_mean_percent",
        compare_yerr_col="recovery_near_true_sem_percent",
    )


def selected_filters(all_filters: Dict[str, float], requested: Optional[Sequence[str]]) -> Dict[str, float]:
    if not requested:
        return dict(all_filters)
    requested_set = {name.upper() for name in requested}
    return {name: rate for name, rate in all_filters.items() if name.upper() in requested_set}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run telescope/filter/exposure sweeps for fit_quality_pipeline.")
    parser.add_argument("--output-root", default=".")
    parser.add_argument("--telescopes", nargs="+", choices=sorted(TELESCOPE_CONFIGS), default=sorted(TELESCOPE_CONFIGS))
    parser.add_argument("--filters", nargs="+", help="Optional filter names to run for every selected telescope.")
    parser.add_argument("--n-samples", type=int, default=100)
    parser.add_argument("--t-exp-min", type=float, default=10.0)
    parser.add_argument("--t-exp-max", type=float, default=1000.0)
    parser.add_argument("--n-exposures", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--nx", type=int, default=fqp.build_default_environment()["nx"])
    parser.add_argument("--ny", type=int, default=fqp.build_default_environment()["ny"])
    parser.add_argument("--example-plots", type=int, default=10, help="Example plots per run, capped at 10.")
    parser.add_argument("--max-retries", type=int, default=7)
    parser.add_argument("--methods", nargs="+", default=["least_squares", "nelder", "leastsq"])
    parser.add_argument("--max-nfev", type=int, default=15000)
    parser.add_argument("--max-workers", type=int, default=12)
    parser.add_argument("--psf-case", help="Use this PSF case from every telescope config.")
    parser.add_argument("--psf-sigma", type=float, help="Override PSF sigma in arcsec for every run.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    output_root = Path(args.output_root).resolve()
    t_exps = exposure_grid(args.t_exp_min, args.t_exp_max, args.n_exposures)

    for telescope in args.telescopes:
        config = load_telescope_config(repo_root, telescope)
        all_filters = valid_sky_filters(config)
        filters_to_run = selected_filters(all_filters, args.filters)
        if not filters_to_run:
            requested = ", ".join(args.filters or [])
            print(f"[skip] {telescope}: no filters with numeric sky background matched {requested!r}")
            continue

        for filter_name, sky_rate in filters_to_run.items():
            print(f"[filter] {telescope}/{filter_name}: {len(t_exps)} exposure point(s)")
            summary_rows: List[Dict[str, Any]] = []
            for t_exp in t_exps:
                run_dir = output_root / telescope / filter_name / exposure_dir_name(float(t_exp))
                env = build_environment(
                    config=config,
                    sky_rate=sky_rate,
                    t_exp=float(t_exp),
                    nx=args.nx,
                    ny=args.ny,
                    psf_case=args.psf_case,
                    psf_sigma_override=args.psf_sigma,
                )
                run_dir.mkdir(parents=True, exist_ok=True)
                write_run_metadata(run_dir, telescope, filter_name, float(t_exp), sky_rate, env)
                print(
                    "[run] "
                    f"{telescope}/{filter_name} t_exp={float(t_exp):.6g}s "
                    f"noise_sigma={env['noise_sigma']:.6g} -> {run_dir}"
                )
                synthetic_manifest, fit_manifest, specs = run_pipeline_point(
                    run_dir=run_dir,
                    env=env,
                    n_samples=args.n_samples,
                    seed=args.seed,
                    example_plots=min(args.example_plots, 10),
                    max_retries=args.max_retries,
                    methods=args.methods,
                    max_nfev=args.max_nfev,
                    max_workers=args.max_workers,
                )
                summary_rows.append(
                    summarize_point(
                        synthetic_manifest=synthetic_manifest,
                        fit_manifest=fit_manifest,
                        specs=specs,
                        telescope=telescope,
                        filter_name=filter_name,
                        t_exp=float(t_exp),
                        sky_rate=sky_rate,
                        noise_sigma=float(env["noise_sigma"]),
                        run_dir=run_dir,
                    )
                )

            results_dir = output_root / telescope / filter_name / "results"
            save_filter_results(pd.DataFrame(summary_rows), results_dir, telescope, filter_name)
            print(f"[results] {telescope}/{filter_name}: {results_dir}")


if __name__ == "__main__":
    main()
