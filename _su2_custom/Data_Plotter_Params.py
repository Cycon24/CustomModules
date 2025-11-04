# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 09:47:28 2025

@author: BriceM
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import re 
import io

import myLogger as MyLog

# ---------------------- User-facing defaults ---------------------- #
DEFAULT_PARAMS_TO_PLOT = [
    "Pressure",
    "Temperature",
    "Mach",
    "Velocity_x",
    "Velocity_y",
    "Velocity_Angle",  # computed here from Velocity_x, Velocity_y
]
LINE_PROBE_FILENAME = "probe_lines.csv"
SURFACE_CSV_FILENAME = "entire_surface_restart.csv"
HISTORY_FILENAME = "history.csv"

base_filepath = "C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\3D_Tests"
Sweep_dir = ["Pb_sweep"]
             # "PARAM_Sweeps_03_6406\\nBlades_sweep",
             # "PARAM_Sweeps_03_6406\\Radius_sweep",
             # "PARAM_Sweeps_02\\NACA_sweep",
             # "PARAM_Sweeps_03_6412\\AoA_sweep",
             # "PARAM_Sweeps_03_6412\\nBlades_sweep",
             # "PARAM_Sweeps_03_6412\\Radius_sweep",]
Image_Type = "PNG"
Image_DPI = 600
# ------------------------------------------------------------------ #


RESIDUALS_HEADERS = ['Time_Iter', 'Outer_Iter', 'Inner_Iter', 'rms[Rho]',
       'rms[RhoU]', 'rms[RhoV]','rms[RhoW]', 'rms[RhoE]', 'rms[nu]']

def _ensure_plots_dir(param_dir: Path) -> Path:
    plots_dir = param_dir / "Plots"
    plots_dir.mkdir(exist_ok=True)
    return plots_dir


def _compute_velocity_angle_deg(vx: np.ndarray, vy: np.ndarray) -> np.ndarray:
    return np.degrees(np.arctan2(vy, vx))


def _read_probe_lines_csv(param_dir: Path) -> pd.DataFrame | None:
    src = param_dir / LINE_PROBE_FILENAME
    if not src.exists():
        print(f"[WARN] Missing {LINE_PROBE_FILENAME} in {param_dir}. Skipping line plots.")
        return None
    return pd.read_csv(src, sep=None, engine="python")


def _read_entire_surface_csv(param_dir: Path) -> pd.DataFrame | None:
    src = param_dir / SURFACE_CSV_FILENAME
    if not src.exists():
        print(f"[WARN] Missing {SURFACE_CSV_FILENAME} in {param_dir}. Skipping contour plots.")
        return None
    return pd.read_csv(src, sep=None, engine="python")


# =============================================================================
# # ---------- history.csv reader + residuals plot ---------- #
# =============================================================================


# --- Robust residuals (history.csv) loader & plotter with smarter iteration detection ---

def _sanitize_colname(name: str) -> str:
    s = str(name).strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in {"'", '"'}:
        s = s[1:-1]
    return s.strip()

def _read_history_csv(param_dir: Path) -> pd.DataFrame | None:
    """
    Read and sanitize history.csv robustly:
    - Normalize newlines and strip trailing commas that cause extra fields
    - Try auto-sep first; fallback to combined regex (comma OR whitespace)
    - Skip malformed lines rather than failing
    - Sanitize headers (strip spaces/quotes)
    """
    src = param_dir / HISTORY_FILENAME
    if not src.exists():
        print(f"[WARN] Missing {HISTORY_FILENAME} in {param_dir}. Skipping residuals plot.")
        return None
    
    df = pd.read_csv(src)
    df.columns = RESIDUALS_HEADERS[:len(df.columns.tolist())]
    

    # # Read raw text & normalize
    # try:
    #     raw = src.read_text(encoding="utf-8-sig", errors="replace")
    # except Exception:
    #     raw = src.read_text(encoding="latin-1", errors="replace")
    # raw = raw.replace("\r\n", "\n").replace("\r", "\n")

    # # Strip trailing commas (phantom empty fields)
    # cleaned_lines = [re.sub(r",\s*$", "", line) for line in raw.split("\n")]
    # cleaned = "\n".join(cleaned_lines)

    # # Try 1: autodetect separator
    # try:
    #     df = pd.read_csv(
    #         io.StringIO(cleaned),
    #         sep=None, engine="python",
    #         skip_blank_lines=True,
    #         on_bad_lines="skip",
    #     )
    # except Exception as e1:
    #     print(f"[INFO] Auto-sep parse failed ({e1}); trying comma-or-whitespace regex...")
    #     # Try 2: comma OR whitespace
    #     try:
    #         df = pd.read_csv(
    #             io.StringIO(cleaned),
    #             sep=r",\s*|\s+",
    #             engine="python",
    #             skip_blank_lines=True,
    #             on_bad_lines="skip",
    #         )
    #     except Exception as e2:
    #         print(f"[WARN] Failed reading {src}: {e2}")
    #         return None

    if df is None or df.empty:
        print(f"[WARN] Empty/invalid {HISTORY_FILENAME} in {param_dir}.")
        return None

    # Clean headers
    # df = df.rename(columns=_sanitize_colname)
    return df

# def _detect_iteration_column(df: pd.DataFrame, user_hint: str | None = None) -> str | None:
#     """
#     Pick the iteration column using (in order):
#       1) user hint (exact match)
#       2) known names (case-insensitive)
#       3) regex patterns
#       4) heuristic: choose the numeric column that looks most like a counter
#          (mostly non-decreasing, many uniques, large range)
#     Returns the ORIGINAL column name (case preserved).
#     """
#     cols = list(df.columns)
#     if not cols:
#         return None

#     # 1) user hint
#     if user_hint is not None:
#         for c in cols:
#             if c == user_hint:
#                 return c
#         # also try case-insensitive match
#         for c in cols:
#             if c.lower() == user_hint.lower():
#                 return c

#     # Normalized map for case-insensitive lookups
#     lower_map = {c.lower(): c for c in cols}

#     # 2) common names (expanded)
#     name_candidates = [
#         "inner_iter", "inneriteration", "inner-iter",
#         "iteration", "iter", "outeriteration", "outer-iter",
#         "time_iter", "timeiteration", "time-iter",
#         "nl_iter", "nonlinear_iter", "nonlineariteration", "nonlinear-iter",
#         "ext_iter", "extiter",
#         "step", "steps", "time_step", "timestep",
#         "it", "i",
#     ]
#     for nc in name_candidates:
#         if nc in lower_map:
#             return lower_map[nc]

#     # 3) regex patterns
#     patterns = [
#         r"^(inner|nonlinear).*(iter|iteration)$",
#         r"^(outer|time).*(iter|iteration)$",
#         r"^iter(ation)?$",
#         r".*(^|_)(time_)?step(s)?$",
#     ]
#     for pat in patterns:
#         regex = re.compile(pat, flags=re.IGNORECASE)
#         for c in cols:
#             if regex.match(c):
#                 return c

#     # 4) heuristic fallback over numeric columns
#     scores = []  # (is_good, frac_nonneg, n_unique, range, colname)
#     for c in cols:
#         s = pd.to_numeric(df[c], errors="coerce")
#         s = s[np.isfinite(s)]
#         if s.size < 5:
#             continue
#         d = np.diff(s)
#         if d.size == 0:
#             continue
#         frac_nonneg = float((d >= 0).mean())
#         n_unique = int(pd.Series(s).nunique())
#         rng = float(np.nanmax(s) - np.nanmin(s))
#         # consider "good" if mostly non-decreasing and reasonably varied
#         is_good = (frac_nonneg >= 0.7) and (n_unique >= 10)
#         scores.append((is_good, frac_nonneg, n_unique, rng, c))

#     if scores:
#         # Prefer "good"; then higher frac_nonneg, then more unique, then larger range
#         scores.sort(key=lambda t: (not t[0], -t[1], -t[2], -t[3]))
#         return scores[0][4]

#     # Give a helpful debug hint
#     print(f"[WARN] Iteration column not found. Available columns: {cols[:12]}{' ...' if len(cols) > 12 else ''}")
#     return None

def _find_rms_columns(df: pd.DataFrame, iter_col: str) -> list[str]:
    # residual cols often look like rms[...]
    cols = []
    for c in df.columns:
        if c == iter_col:
            continue
        cl = c.lower()
        if "rms[" in cl or cl.startswith("rms"):
            cols.append(c)
    return cols

def _is_log10_series(s: pd.Series) -> bool:
    s = pd.to_numeric(s, errors="coerce").dropna()
    return (not s.empty) and ((s < 0).mean() >= 0.2)

def _nice_residual_label(colname: str) -> str:
    m = re.search(r"rms\[(.+?)\]", colname, flags=re.IGNORECASE)
    return m.group(1) if m else colname

def _plot_residuals_history(param_dir: Path, img_format: str, dpi: int, iter_col_hint: str | None = "Inner_Iter"):
    """
    Plot Residuals: RMS vs Iteration from history.csv
    - Iteration on x-axis (auto-picked or user hint)
    - Linear RMS on log y-axis
    """
    df = _read_history_csv(param_dir)
    if df is None:
        return

    # iter_col = _detect_iteration_column(df, user_hint=iter_col_hint)
    # if iter_col is None:
    #     print(f"[WARN] Could not detect iteration column in {param_dir / HISTORY_FILENAME}. Skipping residuals plot.")
    #     return

    # rms_cols = _find_rms_columns(df, iter_col)
    # if not rms_cols:
    #     print(f"[WARN] No residual columns found in {param_dir / HISTORY_FILENAME}. Skipping residuals plot.")
    #     return

    # use = df[[iter_col] + rms_cols].copy()
    # use[iter_col] = pd.to_numeric(use[iter_col], errors="coerce")
    # for c in rms_cols:
    #     use[c] = pd.to_numeric(use[c], errors="coerce")
    # use = use.dropna().sort_values(iter_col, kind="mergesort")
    # if use.empty:
    #     print(f"[WARN] No valid rows in {param_dir / HISTORY_FILENAME} after cleaning.")
    #     return

    plots_dir = _ensure_plots_dir(param_dir)
    fig, ax = plt.subplots()
    
    
    for c in df.columns.tolist():
        if c.startswith('rms'):
            # series = use[c]
            # yvals = np.power(10.0, series.to_numpy(float)) if _is_log10_series(series) else series.to_numpy(float)
            # xvals = use[iter_col].to_numpy(float)
            ax.plot(df["Inner_Iter"], df[c], label=_nice_residual_label(c))

    
    ax.set_title("Residuals")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("RMS")
    # ax.set_yscale("log")  # standard for residuals
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.set_ylim([-9,2])
    ax.minorticks_on()
    ax.legend()

    fig.tight_layout()
    outpath = plots_dir / f"Residuals.{img_format.lower()}"
    fig.savefig(outpath, dpi=dpi, format=img_format.lower())
    plt.close(fig)
    print(f"[OK] {param_dir.name}: saved {outpath.name}")



# =============================================================================
# # --------------------------------------------------------------- #
# =============================================================================

def _validate_columns(df: pd.DataFrame, needed: list[str], context: str) -> bool:
    missing = [c for c in needed if c not in df.columns]
    if missing:
        print(f"[WARN] Missing columns in {context}: {missing}")
        return False
    return True


def _plot_line_property(param_dir: Path, df: pd.DataFrame, prop: str, img_format: str, dpi: int):
    """
    Plot PROPERTY value (x-axis) vs y [in] (y-axis) for each line_id.
    Sort by y (primary) and s_along (secondary). Legend shows line_id and x_ave (inches).
    """
    plots_dir = _ensure_plots_dir(param_dir)

    needed = ["x", "y", "line_id", "s_along", "Velocity_x", "Velocity_y"]
    if not _validate_columns(df, needed, f"{param_dir.name}/{LINE_PROBE_FILENAME}"):
        return

    local_df = df.copy()

    # Compute Velocity_Angle if requested
    if prop == "Velocity_Angle":
        local_df["Velocity_Angle"] = _compute_velocity_angle_deg(
            local_df["Velocity_x"].to_numpy(float),
            local_df["Velocity_y"].to_numpy(float)
        )

    if prop not in local_df.columns:
        print(f"[WARN] Property '{prop}' not found in {LINE_PROBE_FILENAME} for {param_dir.name}. Skipping.")
        return

    fig, ax = plt.subplots()

    for lid in sorted(local_df["line_id"].unique()):
        block = local_df[local_df["line_id"] == lid]
        if block.empty:
            continue

        # Sort by y then s_along for clean connection
        block = block.sort_values(["y", "s_along"], kind="mergesort")

        # PROPERTY on x-axis, y (in) on y-axis
        x_prop = block[prop].to_numpy(float)
        y_in = block["y"].to_numpy(float)

        # Average x-location for this line (already inches in probe CSV)
        x_ave = float(np.nanmean(block["x"].to_numpy(float)))

        ax.plot(x_prop, y_in, label=f"line {lid} (x={x_ave:.1f})")

    ax.set_title(prop)
    ax.set_xlabel(prop if prop != "Velocity_Angle" else "Velocity_Angle [deg]")
    ax.set_ylabel("y [in]")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()
    fig.tight_layout()

    outpath = plots_dir / f"lineplot_{prop}.{img_format.lower()}"
    fig.savefig(outpath, dpi=dpi, format=img_format.lower())
    plt.close(fig)
    print(f"[OK] {param_dir.name}: saved {outpath.name}")


def _plot_contour_property(param_dir: Path, df: pd.DataFrame, prop: str, img_format: str, dpi: int):
    """
    Filled contour of scattered data using triangulation.
    x,y converted from feet to inches. Horizontal colorbar below.
    """
    plots_dir = _ensure_plots_dir(param_dir)

    needed = ["x", "y", "Velocity_x", "Velocity_y"]
    if not _validate_columns(df, needed, f"{param_dir.name}/{SURFACE_CSV_FILENAME}"):
        return

    local_df = df.copy()

    # Coordinates to inches for plotting
    x_in = local_df["x"].to_numpy(float) * 12.0
    y_in = local_df["y"].to_numpy(float) * 12.0

    if prop == "Velocity_Angle":
        local_df["Velocity_Angle"] = _compute_velocity_angle_deg(
            local_df["Velocity_x"].to_numpy(float),
            local_df["Velocity_y"].to_numpy(float)
        )

    if prop not in local_df.columns:
        print(f"[WARN] Property '{prop}' not found in {SURFACE_CSV_FILENAME} for {param_dir.name}. Skipping.")
        return

    z = local_df[prop].to_numpy(float)

    tri = mtri.Triangulation(x_in, y_in)

    fig, ax = plt.subplots()
    tcf = ax.tricontourf(tri, z, levels=50)

    # Horizontal colorbar below
    cbar = fig.colorbar(tcf, ax=ax, orientation="horizontal", pad=0.12)
    cbar.set_label(prop if prop != "Velocity_Angle" else "Velocity_Angle [deg]")

    ax.set_title(prop)
    ax.set_xlabel("x [in]")
    ax.set_ylabel("y [in]")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.set_aspect("equal", adjustable="box")

    fig.tight_layout()
    outpath = plots_dir / f"contour_{prop}.{img_format.lower()}"
    fig.savefig(outpath, dpi=dpi, format=img_format.lower())
    plt.close(fig)
    print(f"[OK] {param_dir.name}: saved {outpath.name}")


def _process_param_folder(param_dir: Path, img_format: str, dpi: int, params_to_plot: list[str]):
    """Create line, contour, and residuals plots for a single Param folder."""
    # Line plots
    probe_df = _read_probe_lines_csv(param_dir)
    if probe_df is not None:
        for prop in params_to_plot:
            _plot_line_property(param_dir, probe_df, prop, img_format, dpi)

    # Contour plots
    surf_df = _read_entire_surface_csv(param_dir)
    if surf_df is not None:
        for prop in params_to_plot:
            _plot_contour_property(param_dir, surf_df, prop, img_format, dpi)

    # NEW: Residuals plot from history.csv
    _plot_residuals_history(param_dir, img_format, dpi)


def _iter_param_dirs(sweep_root: Path):
    for child in sorted(sweep_root.iterdir()):
        if child.is_dir():
            yield child


def _start(sweep_root: str,
    img_format: str = "png",
    dpi: int = 200,
    params_to_plot: list[str] | None = None,
) -> None:
    """
    Actually runs everything
    """
    sr_path = Path(sweep_root).resolve()
    if not sr_path.exists():
        raise FileNotFoundError(f"[ERROR] sweep_root not found: {sr_path}")

    props = DEFAULT_PARAMS_TO_PLOT if params_to_plot is None else list(params_to_plot)

    print(f"[INFO] Sweep root: {sr_path}")
    print(f"[INFO] Output: <Param>/Plots/*.{img_format.lower()} @ {dpi} DPI")
    print(f"[INFO] Properties: {props}")
    print("------------------------------------------------------------")

    any_found = False
    for param_dir in _iter_param_dirs(sr_path):
        any_found = True
        print(f"[INFO] Processing {param_dir.name} ...")
        _process_param_folder(param_dir, img_format=img_format, dpi=dpi, params_to_plot=props)

    if not any_found:
        print(f"[WARN] No subfolders found under {sr_path}")


def main(
    sweep_root: str,
    img_format: str = "png",
    dpi: int = 200,
    params_to_plot: list[str] | None = None,
    log=False
) -> None:
    """
    Generate line, contour, and residual plots for each Param folder in a Sweep.

    Args:
        sweep_root: Path to Sweep directory containing Param subfolders.
        img_format: Image format for saved plots, e.g. 'png', 'pdf', 'svg'. Default 'png'.
        dpi: Figure DPI for saved plots. Default 200.
        params_to_plot: List of property names to plot. If None, uses DEFAULT_PARAMS_TO_PLOT.
        log: If True, wrap in FDTeeLogger and write 'param_plot.log' in sweep_root.
    """
    if log:
        logger = MyLog.FDTeeLogger(filename="param_plot.log", filepath=sweep_root)
        return logger.log(_start, sweep_root, img_format, dpi, params_to_plot)
    else:
        return _start(sweep_root, img_format, dpi, params_to_plot)


if __name__ == "__main__":
    for sweep_dir in Sweep_dir:
        main(base_filepath + "\\" + sweep_dir, Image_Type, Image_DPI, log=False)
