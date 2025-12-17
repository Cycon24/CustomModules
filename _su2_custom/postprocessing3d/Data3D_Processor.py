# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 19:10:11 2025

@author: BriceM
from GPT
"""

# Data3D_Processor.py

from pathlib import Path
from typing import Optional
import warnings

from postprocessing3d.cfd3d_utils import load_main_config
from postprocessing3d.extract_probe_surfaces import extract_probe_surfaces


def _iter_param_dirs(sweep_root: Path):
    """
    Yield immediate subdirectories under sweep_root that look like param folders
    (e.g. AoA_10, AoA_12, ...).
    """
    sweep_root = Path(sweep_root)
    for child in sorted(sweep_root.iterdir()):
        if child.is_dir():
            yield child


def _run_for_param(param_dir: Path, cfg: dict):
    """
    Run probe surface extraction for a single parameter directory.
    Creates:
      - param_dir/Probes/slice_x=..._dx=...csv  (the slabbed cross-sections)
      - param_dir/Logs/post_summary.json/.txt   (diagnostics for each slice)
    """
    print(f"[INFO] 3D Processor: extracting slices for {param_dir.name}")
    extract_probe_surfaces(param_dir, cfg)


def main(
    sweep_root: str,
    cfg_path: str,
    units_path: Optional[str] = None,
) -> None:
    """
    Perform the 3-D preprocessing stage across an entire sweep.

    This is analogous to Data_Processor.main(...) in your 2-D pipeline.

    For each param folder under sweep_root:
      1. reads entire_surface_restart.csv
      2. generates axial "probe surface" slabs at each requested x-station
         (using dx_tolerance_in / min_points_per_slice from the config)
      3. logs slice diagnostics in Logs/post_summary.*

    Args
    ----
    sweep_root : str
        Path to the sweep directory containing subfolders like AoA_10, AoA_12, ...
    cfg_path : str
        Path to post3d_config.yaml
    units_path : Optional[str]
        Currently unused in this stage, but kept for API symmetry with
        Data3D_Plotter_Params.main and future expansion (for example, if we
        someday want unit-aware filtering here).
    """
    sweep_root = Path(sweep_root).resolve()
    if not sweep_root.exists():
        raise FileNotFoundError(f"[ERROR] sweep_root not found: {sweep_root}")

    cfg_path = Path(cfg_path).resolve()
    if not cfg_path.exists():
        raise FileNotFoundError(f"[ERROR] cfg_path not found: {cfg_path}")

    cfg = load_main_config(cfg_path)

    any_found = False
    for param_dir in _iter_param_dirs(sweep_root):
        any_found = True
        _run_for_param(param_dir, cfg)

    if not any_found:
        warnings.warn(f"[WARN] No subfolders found under {sweep_root}")
    else:
        print("[INFO] 3D Processor complete.")


if __name__ == "__main__":
    # Example debug usage
    sweep_example = r"C:\path\to\Sweep\AoA_sweep"
    cfg_example   = r"C:\path\to\post3d_config.yaml"
    main(sweep_example, cfg_example)
