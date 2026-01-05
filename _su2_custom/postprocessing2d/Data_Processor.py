# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 08:42:50 2025

@author: BriceM
"""

# From GPT for sorting through data
"""
Extract near-line probe data from each Param folder under a root "Sweep" directory.

- Reads:   entire_surface_restart.csv  (coords in FEET)
- Lines:   hard-coded below in INCHES
- Tolerance: distance tolerance in INCHES (converted to feet internally)
- Geometry: "ray" (default), "line" (infinite line), or "segment" (bounded)
- Writes:  probe_lines.csv in each Param folder with all original columns + line_id, s_along, dist_to_line
"""

from typing import Iterable
import sys
from pathlib import Path
import warnings
import pandas as pd
import numpy as np

# -----------------------
# USER: define your lines in INCHES here
# Each entry is [x1_in, y1_in, x2_in, y2_in]
# Example:
LINES_INCHES = [
    [3.0, 0.0, 3.0, 1.0],
    [10.0, 0.0, 10.0, 1.0],
    [12.0, 0.0, 12.0, 1.0],
    [15.0, 0.0, 15.0, 1.0],
    [17.0, 0.0, 17.0, 1.0],
    # Add more lines here...
]

Sweep_dir = "PARAMS_RANS_Sweep_01\\NACA_sweep"
p2l_tol_in = 0.1
MIN_LINE_PTS = 25
# -----------------------

EXPECTED_COLUMNS = [
    "PointID","x","y","Density","Momentum_x","Momentum_y","Energy","Pressure",
    "Temperature","Mach","Pressure_Coefficient","Velocity_x","Velocity_y",
    "Laminar_Viscosity","Skin_Friction_Coefficient_x","Skin_Friction_Coefficient_y",
    "Heat_Flux","Y_Plus"
]

def _apply_geometry_t(t_raw: np.ndarray, geometry: str) -> np.ndarray:
    """Apply geometry rule to projection parameter t."""
    if geometry == "segment":
        return np.clip(t_raw, 0.0, 1.0)
    if geometry == "ray":
        return np.maximum(t_raw, 0.0)
    if geometry == "line":
        return t_raw
    raise ValueError(f"Unknown geometry: {geometry!r}. Use 'ray', 'line', or 'segment'.")


def _distances_to_path_all_points(
    x_ft: np.ndarray, y_ft: np.ndarray,
    p1_ft: np.ndarray, v_ft: np.ndarray,
    geometry: str
) -> tuple[np.ndarray, np.ndarray]:
    """
    For ALL points, compute (t_used, d_ft) where:
      t_used = projection parameter after geometry rule,
      d_ft   = perpendicular distance (feet) to the path.
    """
    v2 = v_ft[0]**2 + v_ft[1]**2
    dx = x_ft - p1_ft[0]
    dy = y_ft - p1_ft[1]

    if v2 == 0.0:
        # Degenerate line -> distance to single point p1
        d_ft = np.hypot(dx, dy)
        t_used = np.zeros_like(d_ft)
        return t_used, d_ft

    t_raw = (dx * v_ft[0] + dy * v_ft[1]) / v2
    t_used = _apply_geometry_t(t_raw, geometry)

    px = dx - t_used * v_ft[0]
    py = dy - t_used * v_ft[1]
    d_ft = np.hypot(px, py)
    return t_used, d_ft


def _compute_matches_for_line(
    x_ft: np.ndarray, y_ft: np.ndarray,
    p1_ft: np.ndarray, v_ft: np.ndarray,
    tol_ft: float, geometry: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Vectorized point-to-path filter.
    Returns:
      mask           : bool array of points within tol_ft
      s_along_in     : along-path distance from p1 (inches)
      dist_to_line_in: perpendicular distance (inches)
    """
    t_used, d_ft = _distances_to_path_all_points(x_ft, y_ft, p1_ft, v_ft, geometry)

    mask = d_ft <= tol_ft
    v_mag_in = np.hypot(v_ft[0], v_ft[1]) * 12.0  # feet -> inches
    s_along_in = t_used * v_mag_in
    dist_to_line_in = d_ft * 12.0
    return mask, s_along_in, dist_to_line_in


def _suggest_min_tol_in_for_k(
    x_ft: np.ndarray, y_ft: np.ndarray,
    p1_ft: np.ndarray, v_ft: np.ndarray,
    geometry: str, k: int = 10
) -> float:
    """
    Minimum tolerance (inches) to include at least k points for the given line
    under the selected geometry. This is the k-th smallest distance to path.
    """
    _, d_ft = _distances_to_path_all_points(x_ft, y_ft, p1_ft, v_ft, geometry)
    n = d_ft.size
    if n == 0:
        return float("nan")
    k = max(1, min(k, n))
    kth = np.partition(d_ft, k - 1)[k - 1]
    return float(kth * 12.0)  # inches


def _process_param_folder(
    param_dir: Path,
    lines_inches: Iterable[Iterable[float]],
    tol_in: float,
    geometry: str,
    csv_name: str
) -> None:
    """
    Process one Param folder: read CSV, extract line-probe points, write probe_lines.csv.
    """
    src = param_dir / csv_name
    out = param_dir / "probe_lines.csv"

    if not src.exists():
        warnings.warn(f"[WARN] Missing {csv_name} in {param_dir}. Skipping.")
        return

    # Auto-detect delimiter (comma/tab)
    try:
        df = pd.read_csv(src, sep=None, engine="python")
    except Exception as e:
        print(f"[ERROR] Failed reading {src}: {e}", file=sys.stderr)
        return

    # Validate coord columns
    for col in ("x", "y"):
        if col not in df.columns:
            print(f"[ERROR] Column '{col}' missing in {src}", file=sys.stderr)
            return

    # Keep stable column order where possible
    cols_order = [c for c in EXPECTED_COLUMNS if c in df.columns] + \
                 [c for c in df.columns if c not in EXPECTED_COLUMNS]
    df = df[cols_order]

    # Coordinates from file are FEET
    x_ft = df["x"].to_numpy(dtype=float)
    y_ft = df["y"].to_numpy(dtype=float)

    tol_ft = tol_in / 12.0  # inches -> feet

    all_hits = []
    for line_id, (x1_in, y1_in, x2_in, y2_in) in enumerate(lines_inches):
        # Convert endpoints + direction to FEET
        p1_ft = np.array([x1_in / 12.0, y1_in / 12.0], dtype=float)
        v_ft  = np.array([(x2_in - x1_in) / 12.0, (y2_in - y1_in) / 12.0], dtype=float)

        mask, s_in, d_in = _compute_matches_for_line(
            x_ft, y_ft, p1_ft, v_ft, tol_ft, geometry
        )
        hit_count = int(np.count_nonzero(mask))

        if hit_count < MIN_LINE_PTS:
            min_tol_in = _suggest_min_tol_in_for_k(x_ft, y_ft, p1_ft, v_ft, geometry, k=MIN_LINE_PTS)
            warnings.warn(
                f"[WARN] {param_dir.name}: line_id {line_id} yielded {hit_count} "
                f"points at tol={tol_in:.6g} in; "
                f"min tol for ≥{MIN_LINE_PTS} points ≈ {min_tol_in:.6g} in."
            )

        if hit_count == 0:
            continue

        out_df = df.loc[mask].copy()
        out_df["line_id"] = line_id
        out_df["s_along"] = s_in[mask]          # inches
        out_df["dist_to_line"] = d_in[mask]     # inches
        all_hits.append(out_df)

    # If no hits at all, still write empty with header
    if not all_hits:
        empty_cols = list(df.columns) + ["line_id", "s_along", "dist_to_line"]
        pd.DataFrame(columns=empty_cols).to_csv(out, index=False)
        print(f"[INFO] No points matched lines in {param_dir}. Wrote empty {out.name}.")
        return

    result = pd.concat(all_hits, axis=0, ignore_index=True)

    # Convert output coordinates to INCHES (overwrite x,y)
    if "x" in result.columns:
        result["x"] = result["x"].astype(float) * 12.0
    if "y" in result.columns:
        result["y"] = result["y"].astype(float) * 12.0

    # Sort for usability (by line then along)
    if {"line_id", "s_along"} <= set(result.columns):
        result = result.sort_values(by=["line_id", "s_along"], kind="mergesort").reset_index(drop=True)

    result.to_csv(out, index=False)
    print(f"[OK] Wrote {out}  (rows: {len(result)})")


def _iter_param_dirs(sweep_root: Path):
    """Yield immediate subdirectories (Param folders)."""
    for child in sorted(sweep_root.iterdir()):
        if child.is_dir():
            yield child

def main(
    sweep_root: str,
    lines_inches: list[list[float]] | None = None,
    tol_in: float = 1e-3,
    geometry: str = "ray",
    csv_name: str = "entire_surface_restart.csv",
) -> None:
    """
    Extract near-line probe data (inches output) for each Param folder.

    Args:
        sweep_root: Path to the Sweep directory containing Param subfolders.
        lines_inches: List of [x1_in, y1_in, x2_in, y2_in] lines in INCHES.
                      If None, uses DEFAULT_LINES_INCHES.
        tol_in: Distance tolerance in INCHES (default: 1e-3).
        geometry: 'ray' (default), 'line', or 'segment'.
        csv_name: Input CSV filename inside each Param folder (default: 'entire_surface_restart.csv').
    """
    sr_path = Path(sweep_root).resolve()
    if not sr_path.exists():
        raise FileNotFoundError(f"[ERROR] sweep_root not found: {sr_path}")

    lines = LINES_INCHES if lines_inches is None else list(lines_inches)

    print(f"[INFO] Sweep root: {sr_path}")
    print(f"[INFO] Lines (in): {len(lines)} defined")
    print(f"[INFO] Tolerance: {tol_in} inches   Geometry: {geometry}")
    print("------------------------------------------------------------")

    any_found = False
    for param_dir in _iter_param_dirs(sr_path):
        any_found = True
        _process_param_folder(
            param_dir=param_dir,
            lines_inches=lines,
            tol_in=tol_in,
            geometry=geometry,
            csv_name=csv_name,
        )

    if not any_found:
        warnings.warn(f"[WARN] No subfolders found under {sr_path}")
        
        
if __name__ == "__main__":
    main(Sweep_dir, tol_in=p2l_tol_in)
