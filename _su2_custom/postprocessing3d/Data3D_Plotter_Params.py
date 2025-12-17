# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 19:12:31 2025

@author: BriceM
from GPT
"""

# Data3D_Plotter_Params.py

from pathlib import Path
from typing import Optional
import warnings

from postprocessing3d.cfd3d_utils import (
    load_main_config,
    load_units_config,
)
from postprocessing3d.plot_cross_section_contours import plot_cross_section_contours
from postprocessing3d.plot_radial_averages import plot_radial_averages


def _iter_param_dirs(sweep_root: Path):
    """
    Yield param subdirectories under sweep_root (AoA_10, AoA_12, ...).
    """
    sweep_root = Path(sweep_root)
    for child in sorted(sweep_root.iterdir()):
        if child.is_dir():
            yield child


def _run_for_param(
    param_dir: Path,
    cfg: dict,
    units_cfg: dict,
    img_format: str,
    dpi: int,
):
    """
    For a single parameter folder:
      - Generate contour plots at each axial slice plane
      - Generate radial-averaged line plots
    Output goes in:
      param_dir/Plots/Contours/...
      param_dir/Plots/Radial/...
    """
    print(f"[INFO] 3D Plotter (Params): plotting {param_dir.name}")

    # Contour plots (tricontourf of Mach, Vmag, etc., on (y,z) slice planes)
    plot_cross_section_contours(
        param_dir=param_dir,
        cfg=cfg,
        units_cfg=units_cfg,
    )

    # Radial-averaged line plots (r vs property, one curve per x-station)
    plot_radial_averages(
        param_dir=param_dir,
        cfg=cfg,
        units_cfg=units_cfg,
    )

    # Right now img_format and dpi are not forced down into the plotting
    # helpers. dpi is read from cfg["plotting"]["dpi"], and images save as .png.
    # We keep the arguments here so your external call matches your 2-D style
    # and we can hook these in later without changing your external runner.


def main(
    sweep_root: str,
    cfg_path: str,
    units_path: Optional[str] = None,
    img_format: str = "png",
    dpi: int = 200,
) -> None:
    """
    Produce per-parameter plots across a sweep. Mirrors Data_Plotter_Params.main
    from the 2-D pipeline, but for 3-D data.

    Args
    ----
    sweep_root : str
        Path to the sweep directory with param subfolders (AoA_10, AoA_12, ...)
    cfg_path : str
        Path to post3d_config.yaml
    units_path : Optional[str]
        Optional override path to plot_units.yaml
        (if None, we'll use cfg["plotting"]["units_yaml"])
    img_format : str
        Currently informational; .png is hard-coded in the helpers
    dpi : int
        Currently informational; actual DPI comes from cfg["plotting"]["dpi"]
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

    any_found = False
    for param_dir in _iter_param_dirs(sweep_root):
        any_found = True
        _run_for_param(param_dir, cfg, units_cfg, img_format, dpi)

    if not any_found:
        warnings.warn(f"[WARN] No subfolders found under {sweep_root}")
    else:
        print("[INFO] 3D plotting (Params) complete.")


if __name__ == "__main__":
    # Example debug usage
    sweep_example  = r"C:\path\to\Sweep\AoA_sweep"
    cfg_example    = r"C:\path\to\post3d_config.yaml"
    units_example  = r"C:\path\to\plot_units.yaml"
    main(sweep_example, cfg_example, units_example, img_format="png", dpi=250)
