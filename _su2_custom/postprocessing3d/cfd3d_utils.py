import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
import matplotlib.tri as mtri


# -------------------------------------------------
# Filesystem helpers
# -------------------------------------------------

def ensure_dir(p: Path) -> None:
    """
    Create directory p (and parents) if it doesn't already exist.
    """
    p.mkdir(parents=True, exist_ok=True)


# -------------------------------------------------
# Config helpers
# -------------------------------------------------

def load_main_config(cfg_path: Path) -> dict:
    """
    Load the main post-processing config (post3d_config.yaml).
    """
    import yaml
    return yaml.safe_load(Path(cfg_path).read_text(encoding="utf-8"))


def load_units_config(units_path: Optional[Path], fallback_from_cfg: Optional[Path]) -> dict:
    """
    Load plot_units.yaml-like config. If units_path is None, fall back to the one
    referenced in the main config. If neither exists, return safe defaults.
    """
    import yaml

    if units_path is None:
        if fallback_from_cfg is None:
            return {
                "title": {
                    "template": "{var}: {param} — x={x_in:.2f} in",
                    "template_multi": "{var}: {param}",
                },
                "units": {},
                "limits": {},
                "naming": {
                    "contours_dir": "Plots/Contours",
                    "radial_dir": "Plots/Radial",
                    "probes_dir": "Probes",
                },
            }
        units_path = fallback_from_cfg

    units_path = Path(units_path)
    if not units_path.exists():
        # minimal defaults
        return {
            "title": {
                "template": "{var}: {param} — x={x_in:.2f} in",
                "template_multi": "{var}: {param}",
            },
            "units": {},
            "limits": {},
            "naming": {
                "contours_dir": "Plots/Contours",
                "radial_dir": "Plots/Radial",
                "probes_dir": "Probes",
            },
        }

    return yaml.safe_load(units_path.read_text(encoding="utf-8"))


# -------------------------------------------------
# Reading SU2 surface CSV
# -------------------------------------------------

def read_surface_csv(
    csv_path: Path,
    columns: Dict[str, Any],
) -> pd.DataFrame:
    """
    Read SU2 entire_surface_restart.csv (or similar surface CSV).

    columns mapping example (from YAML):
        x: "x"
        y: "y"
        z: "z"
        density: "Density"
        momentum_x: "Momentum_x"
        momentum_y: "Momentum_y"
        momentum_z: "Momentum_z"
        velocity_x: "Velocity_x"
        velocity_y: "Velocity_y"
        velocity_z: "Velocity_z"
        mach: "Mach"
        pressure: "Pressure"
        temperature: "Temperature"

    Returns DataFrame containing:
        x,y,z,rho,Mx,My,Mz,Vx,Vy,Vz,Mach,Pressure,Temperature
        (some may be derived, see below)
    """
    csv_path = Path(csv_path)

    # Auto-sep read
    df_raw = pd.read_csv(csv_path, sep=None, engine="python")
    df_raw.columns = [c.strip() for c in df_raw.columns]

    def pick(name, default=None):
        return columns.get(name, default)

    # Required coordinate columns
    coord_cols = [pick("x", "x"), pick("y", "y"), pick("z", "z")]
    for c in coord_cols:
        if c not in df_raw.columns:
            raise ValueError(f"Required coordinate column '{c}' not found in {csv_path}")

    # Map into canonical names
    colmap = {
        "x": pick("x", "x"),
        "y": pick("y", "y"),
        "z": pick("z", "z"),
        "rho": pick("density", "Density"),
        "Mx": pick("momentum_x", "Momentum_x"),
        "My": pick("momentum_y", "Momentum_y"),
        "Mz": pick("momentum_z", "Momentum_z"),
        "Vx": pick("velocity_x", "Velocity_x"),
        "Vy": pick("velocity_y", "Velocity_y"),
        "Vz": pick("velocity_z", "Velocity_z"),
        "Mach": pick("mach", "Mach"),
        "Pressure": pick("pressure", "Pressure"),
        "Temperature": pick("temperature", "Temperature"),
    }

    df = pd.DataFrame()
    for out_name, src_name in colmap.items():
        if src_name is None:
            continue
        if src_name in df_raw.columns:
            df[out_name] = df_raw[src_name]

    # derive velocity if missing, using momentum/density
    have_v = all(c in df.columns for c in ["Vx", "Vy", "Vz"])
    have_m = all(c in df.columns for c in ["Mx", "My", "Mz"]) and ("rho" in df.columns)

    if not have_v and have_m:
        rho = df["rho"].replace(0.0, np.nan)
        df["Vx"] = df.get("Vx", df["Mx"] / rho)
        df["Vy"] = df.get("Vy", df["My"] / rho)
        df["Vz"] = df.get("Vz", df["Mz"] / rho)
    elif not have_v and not have_m:
        raise ValueError(
            "Velocity components not present and cannot be derived (missing momentum or density)."
        )
        
   
    # NOTE: exports in ft xyz still
    # Convert to inches
    df["x"] = df["x"].to_numpy() * 12
    df["y"] = df["y"].to_numpy() * 12
    df["z"] = df["z"].to_numpy() * 12
    
    return df


# -------------------------------------------------
# Cylindrical transform about x-axis
# -------------------------------------------------

def add_cylindrical_about_x(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds:
      r        = sqrt(y^2 + z^2)
      theta    = atan2(z, y)
      Vr       = Vy*cos(theta) + Vz*sin(theta)
      Vtheta   = -Vy*sin(theta) + Vz*cos(theta)
      Vmag     = sqrt(Vx^2 + Vy^2 + Vz^2)
      alpha    = atan2(Vtheta, Vx)
      alpha_deg
      Momentum_x, Momentum_r, Momentum_theta (if Mx/My/Mz exist)
    """
    y = df["y"].to_numpy()
    z = df["z"].to_numpy()

    theta = np.arctan2(z, y)
    ct = np.cos(theta)
    st = np.sin(theta)

    Vx = df["Vx"].to_numpy()
    Vy = df["Vy"].to_numpy()
    Vz = df["Vz"].to_numpy()

    Vr = Vy * ct + Vz * st
    Vtheta = -Vy * st + Vz * ct
    Vmag = np.sqrt(Vx * Vx + Vy * Vy + Vz * Vz)
    alpha = np.arctan2(Vtheta, Vx)
    r = np.sqrt(y * y + z * z)

    out = df.copy()
    out["theta"] = theta
    out["r"] = r
    out["Vr"] = Vr
    out["Vtheta"] = Vtheta
    out["Vmag"] = Vmag
    out["alpha"] = alpha
    out["alpha_deg"] = np.degrees(alpha)

    # cylindrical momentum components if momentum is present
    if all(c in out.columns for c in ["Mx", "My", "Mz"]):
        My = out["My"].to_numpy()
        Mz = out["Mz"].to_numpy()
        Mr = My * ct + Mz * st
        Mtheta = -My * st + Mz * ct
        out["Momentum_x"] = out["Mx"]
        out["Momentum_r"] = Mr
        out["Momentum_theta"] = Mtheta

    # Still in ft
    return out


# -------------------------------------------------
# Slab extraction
# -------------------------------------------------

def extract_slab(
    df: pd.DataFrame,
    x0: float,
    dx_tolerance: Optional[float],
    min_points: int,
    max_expand_factor: float = 128.0,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Extract a finite-thickness slab around x0. If dx_tolerance is None,
    adaptively increase dx until we collect at least min_points.

    Returns:
      slice_df : dataframe of points in the slab
      diags    : dict with dx_used, n_points, mean_abs_dx, std_abs_dx, r_min, r_max ...
    """    
    x = df["x"].to_numpy()

    if dx_tolerance is None:
        # start from a tiny fraction of total x-extent
        xr = float(np.nanmax(x) - np.nanmin(x))
        dx_used = max(xr * 1e-4, 1e-6)
        grow = 1.5
        for _ in range(64):
            mask = np.abs(x - x0) <= dx_used / 2.0
            if int(mask.sum()) >= min_points or dx_used > xr * max_expand_factor:
                break
            dx_used *= grow
    else:
        dx_used = float(dx_tolerance)

    mask = np.abs(x - x0) <= dx_used / 2.0
    slice_df = df.loc[mask].copy()

    if len(slice_df) > 0:
        abs_dx = np.abs(slice_df["x"].to_numpy() - x0)
        mean_abs = float(np.mean(abs_dx))
        std_abs = float(np.std(abs_dx, ddof=1)) if len(abs_dx) > 1 else 0.0
        rmin = float(np.nanmin(slice_df["r"])) if "r" in slice_df else float("nan")
        rmax = float(np.nanmax(slice_df["r"])) if "r" in slice_df else float("nan")
    else:
        mean_abs = float("nan")
        std_abs = float("nan")
        rmin = float("nan")
        rmax = float("nan")

    diags = {
        "x0": float(x0),
        "dx_used": float(dx_used),
        "n_points": int(len(slice_df)),
        "mean_abs_dx": mean_abs,
        "std_abs_dx": std_abs,
        "r_min": rmin,
        "r_max": rmax,
    }

    return slice_df, diags


# -------------------------------------------------
# Triangulation for contours
# -------------------------------------------------

def make_triangulation(slice_df: pd.DataFrame) -> mtri.Triangulation:
    """
    Build an unstructured triangulation in the (y,z) plane for contouring.
    """
    y = slice_df["y"].to_numpy()
    z = slice_df["z"].to_numpy()
    return mtri.Triangulation(y, z)


# -------------------------------------------------
# Radial binning stats
# -------------------------------------------------

def radial_bin_stats(
    slice_df: pd.DataFrame,
    r_hub: float,
    r_pipe: float,
    Nr_bins: int,
    variables: List[str],
) -> pd.DataFrame:
    """
    Average each variable on constant-radius bins between r_hub and r_pipe.
    Output columns:
        bin, r_lo, r_hi, r_mid, count, <var>_mean, <var>_std
    """
    r = slice_df["r"].to_numpy()
    bins = np.linspace(r_hub, r_pipe, Nr_bins + 1)
    idx = np.digitize(r, bins) - 1  # 0..Nr_bins-1

    rows = []
    for b in range(Nr_bins):
        sel = (idx == b)
        row = {
            "bin": b,
            "r_lo": bins[b],
            "r_hi": bins[b+1],
            "r_mid": 0.5 * (bins[b] + bins[b+1]),
            "count": int(sel.sum()),
        }
        if np.any(sel):
            chunk = slice_df.loc[sel]
            for q in variables:
                if q not in chunk.columns:
                    row[f"{q}_mean"] = np.nan
                    row[f"{q}_std"] = np.nan
                else:
                    vals = chunk[q].to_numpy()
                    row[f"{q}_mean"] = float(np.nanmean(vals))
                    row[f"{q}_std"] = float(np.nanstd(vals, ddof=1)) if sel.sum() > 1 else 0.0
        else:
            for q in variables:
                row[f"{q}_mean"] = np.nan
                row[f"{q}_std"] = np.nan
        rows.append(row)

    return pd.DataFrame(rows)


# -------------------------------------------------
# Summary logging
# -------------------------------------------------

def append_summary(
    log_json_path: Path,
    log_txt_path: Path,
    param_name: str,
    station_diag: dict,
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Append diagnostics for one slice station to both JSON (cumulative)
    and TXT (line-by-line human readable).
    """
    ensure_dir(log_json_path.parent)

    if log_json_path.exists():
        data = json.loads(log_json_path.read_text(encoding="utf-8") or "{}")
    else:
        data = {}

    entry = dict(station_diag)
    if extra:
        entry.update(extra)

    data.setdefault(param_name, [])
    data[param_name].append(entry)

    log_json_path.write_text(json.dumps(data, indent=2), encoding="utf-8")

    with open(log_txt_path, "a", encoding="utf-8") as f:
        f.write(
            f"[{param_name}] "
            f"x0={entry.get('x0'):.6g}, "
            f"dx={entry.get('dx_used'):.6g}, "
            f"n={entry.get('n_points')}, "
            f"mean|x-x0|={entry.get('mean_abs_dx'):.6g}, "
            f"std={entry.get('std_abs_dx'):.6g}, "
            f"r_min={entry.get('r_min'):.6g}, "
            f"r_max={entry.get('r_max'):.6g}\n"
        )


# -------------------------------------------------
# Titles and names
# -------------------------------------------------

def title_from_template(template: str, var: str, param: str, x_in: float) -> str:
    """
    Safely format plot titles with {var}, {param}, {x_in}.
    """
    try:
        return template.format(var=var, param=param, x_in=x_in)
    except Exception:
        return f"{var}: {param} — x={x_in:.2f} in"


def infer_param_name(path_like: Path) -> str:
    """
    Infer param (sweep condition) name from folder name, e.g. AoA_10.
    """
    return Path(path_like).name
