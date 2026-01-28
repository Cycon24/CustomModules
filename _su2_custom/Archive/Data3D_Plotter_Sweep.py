# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 19:13:02 2025

@author: BriceM
"""

# Data3D_Plotter_Sweep.py

from pathlib import Path
from typing import Optional, Dict, List
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    from cfd3d_utils import (
        load_main_config,
        load_units_config,
        radial_bin_stats,
        ensure_dir,
    )
else:
    from postprocessing3d.cfd3d_utils import (
        load_main_config,
        load_units_config,
        radial_bin_stats,
        ensure_dir,
    )


def _iter_param_dirs(sweep_root: Path):
    """
    Yield param subdirectories under sweep_root (e.g., AoA_10, AoA_12, ...).
    """
    sweep_root = Path(sweep_root)
    for child in sorted(sweep_root.iterdir()):
        if child.is_dir():
            yield child


def _load_slice_for_station(param_dir: Path, x_target: float) -> pd.DataFrame:
    """
    Find and load the probe slice CSV from param_dir/Probes whose filename
    encodes slice_x=<x_target>. Returns None if not found.

    We look for files named like slice_x=12.000000_dx=0.050000.csv.
    We'll pick the one with x value closest to x_target.
    """
    probe_dir = Path(param_dir) / "Probes"
    if not probe_dir.exists():
        return None

    best_df = None
    best_err = None
    best_name = None

    for f in probe_dir.glob("slice_x=*_dx=*.csv"):
        try:
            x_in = float(f.stem.split("slice_x=")[1].split("_dx=")[0])
        except Exception:
            continue
        err = abs(x_in - x_target)
        if (best_err is None) or (err < best_err):
            try:
                cand_df = pd.read_csv(f)
            except Exception:
                continue
            best_err = err
            best_df = cand_df
            best_name = f.name

    return best_df


def _radial_stats_for_slice(df_slice: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Run radial_bin_stats on a single slice df using global geometry/bin settings.
    Returns DataFrame with columns:
      bin,r_lo,r_hi,r_mid,count,<var>_mean,<var>_std,...
    """
    r_hub = float(cfg["geometry"]["r_hub"])
    r_pipe = float(cfg["geometry"]["r_pipe"])
    Nr = int(cfg["radial_averages"]["Nr_bins"])
    vars_for_radial = list(cfg["radial_averages"]["variables"])

    return radial_bin_stats(df_slice, r_hub, r_pipe, Nr, vars_for_radial)


def _plot_sweep_overlay_for_station(
    sweep_root: Path,
    cfg: dict,
    units_cfg: dict,
    x_target: float,
    img_format: str,
    dpi: int,
):
    """
    For a given axial station x_target:
      - for each param folder in sweep_root
      - load its slice closest to x_target
      - compute radial bin stats
      - overlay curves from each param for each variable

    Save:
      SweepPlots/x=<x_target>/<var>_sweepCompare.png
    Also dump a CSV with stacked data across params.
    """
    sweep_root = Path(sweep_root)
    out_root = sweep_root / "SweepPlots"
    ensure_dir(out_root)

    # get param dirs
    param_dirs = list(_iter_param_dirs(sweep_root))
    if not param_dirs:
        warnings.warn(f"[WARN] No param subfolders under {sweep_root} for sweep plotting.")
        return

    # build {param_name: radial_stats_df} for this x_target
    per_param_stats = {}
    for pdir in param_dirs:
        param_name = pdir.name
        df_slice = _load_slice_for_station(pdir, x_target)
        if df_slice is None or len(df_slice) == 0:
            warnings.warn(f"[WARN] {param_name}: no slice data near x={x_target}")
            continue
        stats = _radial_stats_for_slice(df_slice, cfg)
        per_param_stats[param_name] = stats

    if not per_param_stats:
        warnings.warn(f"[WARN] No usable slice data at x={x_target} across sweep.")
        return

    # everyone should share the same r_mid; just grab one
    sample_param = next(iter(per_param_stats))
    r_mid = per_param_stats[sample_param]["r_mid"].to_numpy()

    # variables we will overlay (same as radial_averages.variables)
    vars_for_radial = list(cfg["radial_averages"]["variables"])

    # output dir for this station
    station_dir = out_root / f"x={x_target:.2f}"
    ensure_dir(station_dir)

    for var in vars_for_radial:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)

        # CSV stack builder for this var
        stack_rows = []

        for param_name, stats_df in per_param_stats.items():
            if f"{var}_mean" not in stats_df.columns:
                continue

            yvals = stats_df[f"{var}_mean"].to_numpy()
            ax.plot(r_mid, yvals, label=param_name)

            # accumulate rows for CSV
            tmp = stats_df[["r_mid", f"{var}_mean", f"{var}_std", "count"]].copy()
            tmp["param"] = param_name
            tmp["x_target"] = x_target
            stack_rows.append(tmp)

        unit_lbl = units_cfg.get("units", {}).get(var, "")
        ax.set_xlabel("r [in]")
        ax.set_ylabel(f"{var} {'[' + unit_lbl + ']' if unit_lbl else ''}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

        # Title for sweep comparison:
        # we'll reuse template_multi but add the station to make it obvious
        tmpl_multi = units_cfg.get("title", {}).get(
            "template_multi",
            "{var}: {param}",
        )
        # For sweep overlay, {param} doesn't quite make sense (there are many),
        # so we'll fill it with something descriptive.
        sweep_label = f"AllParams @ x={x_target:.2f} in"
        title = tmpl_multi.format(var=var, param=sweep_label)
        ax.set_title(title)

        out_img = station_dir / f"{var}_sweepCompare.{img_format}"
        fig.savefig(out_img, bbox_inches="tight")
        plt.close(fig)

        # Save CSV with multi-param comparison for this var @ this station
        if stack_rows:
            stack_df = pd.concat(stack_rows, ignore_index=True)
            out_csv = station_dir / f"{var}_sweepCompare.csv"
            stack_df.to_csv(out_csv, index=False)


def main(
    sweep_root: str,
    cfg_path: str,
    units_path: Optional[str] = None,
    img_format: str = "png",
    dpi: int = 200,
) -> None:
    """
    Sweep-level comparison across different parameter cases.

    This is the 3-D analog of Data_Plotter_Sweep.main for your 2-D pipeline.

    For each axial slice location (each x-station in cfg["slices"]["x_stations_in"]):
      - load each param folder's nearest slice_x=... file
      - compute radial bin averages
      - overlay all params on one plot for each variable in
        cfg["radial_averages"]["variables"]

    Plots & CSVs are saved under:
      <sweep_root>/SweepPlots/x=<station>/

    Args
    ----
    sweep_root : str
        Path to the sweep directory (contains AoA_10/, AoA_12/, ...)
    cfg_path : str
        Path to post3d_config.yaml
    units_path : Optional[str]
        Optional override to plot_units.yaml
    img_format : str
        Image format extension to save, e.g. "png"
    dpi : int
        DPI for figures
    """
    sweep_root = Path(sweep_root).resolve()
    if not sweep_root.exists():
        raise FileNotFoundError(f"[ERROR] sweep_root not found: {sweep_root}")

    cfg_path = Path(cfg_path).resolve()
    if not cfg_path.exists():
        raise FileNotFoundError(f"[ERROR] cfg_path not found: {cfg_path}")

    cfg = load_main_config(cfg_path)
    units_cfg = load_units_config(
        units_path=Path(units_path) if units_path else None,
        fallback_from_cfg=Path(cfg["plotting"]["units_yaml"]) if "plotting" in cfg else None,
    )

    # stations we care about come straight from config
    stations = list(cfg["slices"]["x_stations_in"])

    if not stations:
        warnings.warn("[WARN] No x_stations_in defined in config for sweep plotting.")
        return

    print("[INFO] 3D Plotter (Sweep): generating sweep comparison plots...")
    for x_target in stations:
        _plot_sweep_overlay_for_station(
            sweep_root=sweep_root,
            cfg=cfg,
            units_cfg=units_cfg,
            x_target=float(x_target),
            img_format=img_format,
            dpi=dpi,
        )

    print("[INFO] 3D plotting (Sweep) complete.")


if __name__ == "__main__":
    # Example debug usage
    sweep_example  = r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\AoA_rt_sweep\AoA_rt_2_4"
    cfg_example    = 'post3d_config.yaml'
    units_example  = r'C:\Users\BriceM\Documents\Modules\_su2_custom\3dpostprocessing\plot_units.yaml'
    main(sweep_example, cfg_example, units_example, img_format="png", dpi=250)
