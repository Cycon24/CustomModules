# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 10:31:59 2025

@author: BriceM
"""

# Again from GPT
"""
sweep_line_comparisons.py  — with CSV export

Callable entry point:
    main(sweep_root: str, img_format: str = "png", dpi: int = 200,
         params_to_plot: list[str] | None = None) -> None

Creates sweep-level comparison plots AND CSVs:
- For each line_id and for each property:
  • Plot PROPERTY (x-axis) vs y [in] (y-axis), one curve per Param folder.
  • Save a tidy CSV with the sorted data that feeds the plot.

Outputs:
- Figures: <sweep_root>/Plots/line_{line_id}_x_{x_ave:.1f}/comparison_{Property}.<ext>
- CSV:     <sweep_root>/Plots/line_{line_id}_x_{x_ave:.1f}/comparison_{Property}.csv

Assumptions:
- Each Param folder under sweep_root contains a probe_lines.csv with columns:
  x [in], y [in], line_id, s_along [in], dist_to_line [in], Velocity_x, Velocity_y,
  plus properties like Pressure, Temperature, Mach, etc.
- Velocity_Angle is computed here as atan2(Vy, Vx) in degrees.
"""

from pathlib import Path
import warnings
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import yaml  # PyYAML
except Exception:
    yaml = None  # Handle gracefully if not installed



PROBE_FILE = "probe_lines.csv"

DEFAULT_PARAMS_TO_PLOT = [
    "Pressure",
    "Temperature",
    "Mach",
    "Velocity_x",
    "Velocity_y",
    "Velocity_Angle",  # computed from Velocity_x, Velocity_y (degrees)
]

Sweep_dir = "PARAM_Sweeps_02\\NACA_sweep"
Image_Type = "PNG"
Image_DPI = 200


def _load_subs_yaml(path: Optional[str]) -> dict:
    if path is None:
        return {}
    if yaml is None:
        raise RuntimeError("PyYAML is required for substitutions. Install with: pip install pyyaml")
    cf = Path(path)
    if not cf.exists():
        raise FileNotFoundError(f"Substitution YAML not found: {cf}")
    with cf.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if "params" not in data:
        data = {"params": data}
    if "defaults" not in data:
        data["defaults"] = {}
    data["defaults"].setdefault("wrap_math", True)
    data["defaults"].setdefault("string_mode", "mathrm")
    data["defaults"].setdefault("joiner", ";\\ ")
    return data

def _escape_for_mathmathrm(s: str) -> str:
    return s.replace("_", r"\_")

def _try_float(val: str):
    try:
        return True, float(val)
    except Exception:
        return False, 0.0

def _format_component(param: str, value: str, cfg: dict, defaults: dict) -> Tuple[str, bool]:
    symbol = cfg.get("symbol", param)
    wrap_math = cfg.get("wrap_math", defaults.get("wrap_math", True))
    numeric_cfg = cfg.get("numeric", {})
    string_cfg = cfg.get("string", {})
    string_mode = defaults.get("string_mode", "mathrm")

    is_num, fval = _try_float(value)
    if is_num:
        fmt = numeric_cfg.get("fmt", "")
        try:
            value_str = format(fval, fmt) if fmt else str(fval)
        except Exception:
            value_str = str(fval)
        units = numeric_cfg.get("units", "")
        comp = f"{symbol}={value_str}{units}"
    else:
        mapped = string_cfg.get("map", {}).get(value, None)
        if mapped is not None:
            comp = f"{symbol}={mapped}"
        else:
            if string_mode == "mathrm":
                comp = f"{symbol}=\\mathrm{{{_escape_for_mathmathrm(value)}}}"
            else:
                comp = f"{symbol}={value}"
    return comp, wrap_math

def _parse_pairs(folder_name: str) -> List[Tuple[str, str]]:
    toks = folder_name.split("_")
    pairs = []
    for i in range(0, len(toks) - 1, 2):
        pairs.append((toks[i], toks[i+1]))
    return pairs

def _make_pretty_label(folder_name: str, subs: dict) -> str:
    pairs = _parse_pairs(folder_name)
    if not pairs:
        return folder_name

    params_cfg = subs.get("params", {})
    defaults = subs.get("defaults", {"wrap_math": True, "string_mode": "mathrm", "joiner": ";\\ "})

    parts, any_math = [], False
    for p, v in pairs:
        cfg = params_cfg.get(p, {})
        comp, wrap = _format_component(p, v, cfg, defaults)
        parts.append(comp)
        any_math = any_math or wrap

    joined = (defaults.get("joiner", ";\\ ")).join(parts)
    return f"${joined}$" if any_math else joined

# ---------------- Core data utilities (unchanged logic) ---------------- #

def _read_probe_df(param_dir: Path) -> pd.DataFrame | None:
    src = param_dir / PROBE_FILE
    if not src.exists():
        warnings.warn(f"[WARN] Missing {PROBE_FILE} in {param_dir.name}; skipping this folder.")
        return None
    try:
        df = pd.read_csv(src, sep=None, engine="python")
    except Exception as e:
        warnings.warn(f"[WARN] Failed to read {src}: {e}")
        return None
    needed = {"x", "y", "line_id", "s_along", "Velocity_x", "Velocity_y"}
    missing = [c for c in needed if c not in df.columns]
    if missing:
        warnings.warn(f"[WARN] {param_dir.name}/{PROBE_FILE} missing columns: {missing}")
        return None
    return df

def _compute_velocity_angle_deg(vx: np.ndarray, vy: np.ndarray) -> np.ndarray:
    return np.degrees(np.arctan2(vy, vx))

def _label_for_prop(prop: str) -> str:
    return "Velocity_Angle [deg]" if prop == "Velocity_Angle" else prop

def _iter_param_dirs(sweep_root: Path):
    for child in sorted(sweep_root.iterdir()):
        if child.is_dir():
            yield child

def _ensure_line_folder(root_plots: Path, line_id: int, x_ave_in: float) -> Path:
    subdir = root_plots / f"line_{line_id}_x_{x_ave_in:.1f}"
    subdir.mkdir(parents=True, exist_ok=True)
    return subdir

def _collect_data(sweep_root: Path) -> Dict[str, pd.DataFrame]:
    data: Dict[str, pd.DataFrame] = {}
    for param_dir in _iter_param_dirs(sweep_root):
        df = _read_probe_df(param_dir)
        if df is not None and not df.empty:
            data[param_dir.name] = df
    if not data:
        warnings.warn(f"[WARN] No valid {PROBE_FILE} found under {sweep_root}")
    return data

def _common_line_ids(data: Dict[str, pd.DataFrame]) -> List[int]:
    sets = []
    for _, df in data.items():
        lids = set(int(x) for x in df["line_id"].unique())
        sets.append(lids)
    if not sets:
        return []
    intersection = set.intersection(*sets)
    union = set.union(*sets)
    if union != intersection:
        warnings.warn(
            f"[WARN] Not all folders share identical line_ids. "
            f"Using intersection only; missing in some folders: {sorted(union - intersection)}"
        )
    return sorted(intersection)

def _global_x_ave_for_line(data: Dict[str, pd.DataFrame], line_id: int) -> float:
    vals = []
    for df in data.values():
        blk = df[df["line_id"] == line_id]
        if not blk.empty:
            vals.append(float(np.nanmean(blk["x"].to_numpy(float))))
    if not vals:
        return float("nan")
    x_mean = float(np.nanmean(vals))
    x_min, x_max = float(np.nanmin(vals)), float(np.nanmax(vals))
    if (x_max - x_min) > 0.1:
        warnings.warn(
            f"[WARN] x_ave spread for line {line_id} across folders is {x_max - x_min:.3f} in "
            f"(min={x_min:.3f}, max={x_max:.3f}). Using mean {x_mean:.3f} in naming."
        )
    return x_mean

def _build_tidy_dataframe_for_line_and_prop(data: Dict[str, pd.DataFrame], line_id: int, prop: str) -> pd.DataFrame:
    rows = []
    for folder_name, df in sorted(data.items()):
        blk = df[df["line_id"] == line_id]
        if blk.empty:
            continue
        if prop == "Velocity_Angle":
            prop_vals = _compute_velocity_angle_deg(
                blk["Velocity_x"].to_numpy(float),
                blk["Velocity_y"].to_numpy(float)
            )
        else:
            if prop not in blk.columns:
                warnings.warn(f"[WARN] {folder_name} missing '{prop}'. Skipping.")
                continue
            prop_vals = blk[prop].to_numpy(float)
        blk = blk.assign(_prop=prop_vals).sort_values(["y", "s_along"], kind="mergesort")
        rows.append(pd.DataFrame({
            "folder": folder_name,
            "line_id": int(line_id),
            "y_in": blk["y"].to_numpy(float),
            "s_along_in": blk["s_along"].to_numpy(float),
            "x_in": blk["x"].to_numpy(float),
            prop: blk["_prop"].to_numpy(float),
        }))
    if not rows:
        return pd.DataFrame(columns=["folder","line_id","y_in","s_along_in","x_in",prop])
    tidy = pd.concat(rows, axis=0, ignore_index=True)
    return tidy  # keep raw; plotting will order the groups

# ---------------- New: numeric sort helpers ---------------- #

def _extract_sort_value(folder_name: str, sort_key_param: Optional[str] = None) -> Tuple[int, float, str]:
    """
    Returns a tuple sort key:
      (is_non_numeric, numeric_value, folder_name)
    - is_non_numeric: 0 if numeric, 1 otherwise (so numerics come first).
    - numeric_value:  float value if numeric, else +inf (not used).
    - folder_name:    used as a stable tiebreaker / alpha order.
    Strategy:
      * If sort_key_param given: use its value if present & numeric.
      * Else: try first pair; if non-numeric, use the first numeric among pairs.
      * If nothing numeric -> mark as non-numeric.
    """
    pairs = _parse_pairs(folder_name)
    if not pairs:
        return (1, float("inf"), folder_name)

    def val_for_param(pairs, key):
        for p, v in pairs:
            if p == key:
                ok, f = _try_float(v)
                if ok: return f
                return None
        return None

    # 1) explicit key
    if sort_key_param:
        f = val_for_param(pairs, sort_key_param)
        if f is not None:
            return (0, f, folder_name)

    # 2) first pair
    ok, f = _try_float(pairs[0][1])
    if ok:
        return (0, f, folder_name)

    # 3) first numeric among pairs
    for _, v in pairs[1:]:
        ok, f = _try_float(v)
        if ok:
            return (0, f, folder_name)

    # 4) no numeric
    return (1, float("inf"), folder_name)

def _ordered_folders_for_plot(tidy: pd.DataFrame, sort_by_value: bool, sort_key_param: Optional[str], ascending: bool) -> List[str]:
    folders = list(tidy["folder"].unique())
    if not sort_by_value:
        return sorted(folders)  # default alphabetical

    keys = [(_extract_sort_value(name, sort_key_param), name) for name in folders]
    # keys: [((is_non_numeric, numeric_value, folder_name), name), ...]
    keys.sort(key=lambda t: (t[0][0], t[0][1], t[0][2]))  # numerics first, asc by value, then name
    ordered = [name for _, name in keys]
    if not ascending:
        # Keep non-numerics at the end, but reverse the numeric block
        numeric_names = [n for ((nn, _, _), n) in keys if nn == 0]
        non_numeric_names = [n for ((nn, _, _), n) in keys if nn == 1]
        return list(reversed(numeric_names)) + non_numeric_names
    return ordered

# ---------------- Plotting ---------------- #

def _plot_comparison_for_line_and_prop(
    sweep_root: Path,
    root_plots: Path,
    tidy: pd.DataFrame,
    line_id: int,
    prop: str,
    x_ave_in: float,
    img_format: str,
    dpi: int,
    pretty_labels: Dict[str, str],
    sort_by_value: bool,
    sort_key_param: Optional[str],
    ascending: bool,
):
    line_folder = _ensure_line_folder(root_plots, line_id, x_ave_in if np.isfinite(x_ave_in) else 0.0)
    fig, ax = plt.subplots()

    # decide plotting/legend order
    folder_order = _ordered_folders_for_plot(tidy, sort_by_value, sort_key_param, ascending)
    n_series = len(folder_order)
    legend_kwargs = dict(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    legend_ncol = 2 if n_series > 10 else 1

    have_any = False
    for folder_name in folder_order:
        blk = tidy[tidy["folder"] == folder_name]
        if blk.empty:
            continue
        x_vals = blk[prop].to_numpy(float)
        y_vals = blk["y_in"].to_numpy(float)
        label = pretty_labels.get(folder_name, folder_name)
        ax.plot(x_vals, y_vals, label=label)
        have_any = True

    if not have_any:
        plt.close(fig)
        return

    title = f"Line {line_id} (x={x_ave_in:.1f}): {prop}"
    ax.set_title(title)
    ax.set_xlabel("Velocity_Angle [deg]" if prop == "Velocity_Angle" else prop)
    ax.set_ylabel("y [in]")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend(ncol=legend_ncol, **legend_kwargs)
    fig.subplots_adjust(right=0.75)

    outpath = line_folder / f"comparison_{prop}.{img_format.lower()}"
    fig.savefig(outpath, dpi=dpi, format=img_format.lower(), bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved {outpath.relative_to(sweep_root)}")

# ---------------- Public API ---------------- #

def main(
    sweep_root: str,
    img_format: str = "png",
    dpi: int = 200,
    params_to_plot: List[str] | None = None,
    subs_yaml: Optional[str] = None,
    sort_by_value: bool = True,
    sort_key_param: Optional[str] = None,
    ascending: bool = True,
) -> None:
    """
    Create comparison plots across Param folders for each (line_id, property),
    with optional legend substitutions (YAML) and numeric ordering.

    Args:
        sweep_root: path to Sweep directory
        img_format: 'png' | 'pdf' | 'svg' | ...
        dpi: integer DPI
        params_to_plot: properties list or None for default
        subs_yaml: path to YAML substitutions (optional)
        sort_by_value: if True, plot/legend order is by numeric value (low→high)
        sort_key_param: which Param to sort by; if None, auto-select as described
        ascending: True (low→high) or False (high→low)
    """
    sr_path = Path(sweep_root).resolve()
    if not sr_path.exists():
        raise FileNotFoundError(f"[ERROR] sweep_root not found: {sr_path}")

    props = DEFAULT_PARAMS_TO_PLOT if params_to_plot is None else list(params_to_plot)

    data = _collect_data(sr_path)
    if not data:
        return

    line_ids = _common_line_ids(data)
    if not line_ids:
        warnings.warn("[WARN] No common line_ids found; nothing to plot.")
        return

    subs = _load_subs_yaml(subs_yaml) if subs_yaml else {}
    pretty_labels = {k: _make_pretty_label(k, subs) for k in data.keys()} if subs else {k: k for k in data.keys()}

    root_plots = sr_path / "Plots"
    root_plots.mkdir(exist_ok=True)

    print(f"[INFO] Sweep root: {sr_path}")
    print(f"[INFO] Folders: {len(data)}  |  line_ids: {line_ids}")
    print(f"[INFO] Using substitutions: {'yes' if subs else 'no'}")
    print(f"[INFO] Sorting by value: {sort_by_value}  (key={sort_key_param or 'auto'}, asc={ascending})")
    print(f"[INFO] Saving into: {root_plots}")
    print(f"[INFO] Format: .{img_format.lower()}  |  DPI: {dpi}")
    print("------------------------------------------------------------")

    for lid in line_ids:
        x_ave_in = _global_x_ave_for_line(data, lid)
        for prop in props:
            tidy = _build_tidy_dataframe_for_line_and_prop(data, lid, prop)
            if tidy.empty:
                continue
            _plot_comparison_for_line_and_prop(
                sweep_root=sr_path,
                root_plots=root_plots,
                tidy=tidy,
                line_id=lid,
                prop=prop,
                x_ave_in=x_ave_in,
                img_format=img_format,
                dpi=dpi,
                pretty_labels=pretty_labels,
                sort_by_value=sort_by_value,
                sort_key_param=sort_key_param,
                ascending=ascending,
            )



if __name__ == "__main__":
    main(Sweep_dir, Image_Type, Image_DPI, subs_yaml="label_subs.yaml")