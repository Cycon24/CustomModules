import os
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
import matplotlib.tri as mtri

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def read_surface_csv(
    csv_path: Path,
    columns: Dict[str, Any],
    usecols: Optional[List[str]] = None,
    chunksize: Optional[int] = None,
) -> pd.DataFrame:
    csv_path = Path(csv_path)
    if usecols is None:
        with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline().strip()
    df = pd.read_csv(csv_path, sep=None, engine="python")
    df.columns = [c.strip() for c in df.columns]

    def pick(name, default=None):
        return columns.get(name, default)

    req = [pick("x", "x"), pick("y", "y"), pick("z", "z")]
    for c in req:
        if c not in df.columns:
            raise ValueError(f"Required coordinate column '{c}' not found in {csv_path}")

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

    out = pd.DataFrame()
    for k, src in colmap.items():
        if src is None:
            continue
        if src in df.columns:
            out[k] = df[src]

    have_v = all(x in out.columns for x in ["Vx","Vy","Vz"])
    have_m = all(x in out.columns for x in ["Mx","My","Mz"]) and ("rho" in out.columns)

    if not have_v and have_m:
        rho = out["rho"].replace(0.0, np.nan)
        out["Vx"] = out.get("Vx", out["Mx"] / rho)
        out["Vy"] = out.get("Vy", out["My"] / rho)
        out["Vz"] = out.get("Vz", out["Mz"] / rho)
    elif not have_v and not have_m:
        raise ValueError("Velocity components not present and cannot be derived (missing momentum or density).")

    for extra in ["Mx","My","Mz","Mach","Pressure","Temperature"]:
        if extra not in out.columns and colmap.get(extra) in df.columns:
            out[extra] = df[colmap[extra]]

    return out

def add_cylindrical_about_x(df: pd.DataFrame) -> pd.DataFrame:
    y = df["y"].to_numpy()
    z = df["z"].to_numpy()
    theta = np.arctan2(z, y)
    ct = np.cos(theta)
    st = np.sin(theta)

    Vx = df["Vx"].to_numpy()
    Vy = df["Vy"].to_numpy()
    Vz = df["Vz"].to_numpy()

    Vr = Vy*ct + Vz*st
    Vtheta = -Vy*st + Vz*ct
    Vmag = np.sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
    alpha = np.arctan2(Vtheta, Vx)
    r = np.sqrt(y*y + z*z)

    df = df.copy()
    df["theta"] = theta
    df["r"] = r
    df["Vr"] = Vr
    df["Vtheta"] = Vtheta
    df["Vmag"] = Vmag
    df["alpha"] = alpha
    df["alpha_deg"] = np.degrees(alpha)

    if all(c in df.columns for c in ["Mx","My","Mz"]):
        My = df["My"].to_numpy()
        Mz = df["Mz"].to_numpy()
        Mr = My*ct + Mz*st
        Mtheta = -My*st + Mz*ct
        df["Momentum_x"] = df["Mx"]
        df["Momentum_r"] = Mr
        df["Momentum_theta"] = Mtheta

    return df

def extract_slab(
    df: pd.DataFrame,
    x0: float,
    dx_tolerance: Optional[float],
    min_points: int,
    max_expand_factor: float = 128.0,
):
    x = df["x"].to_numpy()
    if dx_tolerance is None:
        xr = float(np.nanmax(x) - np.nanmin(x))
        dx = max(xr * 1e-4, 1e-6)
        factor = 1.5
        iters = 0
        while True:
            mask = np.abs(x - x0) <= dx/2.0
            n = int(mask.sum())
            if n >= min_points or dx > xr * max_expand_factor:
                break
            dx *= factor
            iters += 1
            if iters > 64:
                break
    else:
        dx = float(dx_tolerance)

    mask = np.abs(x - x0) <= dx/2.0
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
        "dx_used": float(dx),
        "n_points": int(len(slice_df)),
        "mean_abs_dx": mean_abs,
        "std_abs_dx": std_abs,
        "r_min": rmin,
        "r_max": rmax,
    }
    return slice_df, diags

def make_triangulation(slice_df: pd.DataFrame):
    y = slice_df["y"].to_numpy()
    z = slice_df["z"].to_numpy()
    tri = mtri.Triangulation(y, z)
    return tri

def radial_bin_stats(
    slice_df: pd.DataFrame,
    r_hub: float,
    r_pipe: float,
    Nr_bins: int,
    variables: List[str],
) -> pd.DataFrame:
    r = slice_df["r"].to_numpy()
    bins = np.linspace(r_hub, r_pipe, Nr_bins+1)
    idx = np.digitize(r, bins) - 1

    rows = []
    for b in range(Nr_bins):
        sel = (idx == b)
        if not np.any(sel):
            row = {"bin": b, "r_lo": bins[b], "r_hi": bins[b+1], "r_mid": 0.5*(bins[b]+bins[b+1]), "count": 0}
            for q in variables:
                row[f"{q}_mean"] = np.nan
                row[f"{q}_std"] = np.nan
            rows.append(row)
            continue

        chunk = slice_df.loc[sel]
        row = {"bin": b, "r_lo": bins[b], "r_hi": bins[b+1], "r_mid": 0.5*(bins[b]+bins[b+1]), "count": int(sel.sum())}
        for q in variables:
            if q not in chunk.columns:
                row[f"{q}_mean"] = np.nan
                row[f"{q}_std"] = np.nan
            else:
                vals = chunk[q].to_numpy()
                row[f"{q}_mean"] = float(np.nanmean(vals))
                row[f"{q}_std"] = float(np.nanstd(vals, ddof=1)) if sel.sum()>1 else 0.0
        rows.append(row)

    return pd.DataFrame(rows)

def append_summary(
    log_json_path: Path,
    log_txt_path: Path,
    param_name: str,
    station_diag: dict,
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    ensure_dir(log_json_path.parent)
    if log_json_path.exists():
        import json as _json
        data = _json.loads(log_json_path.read_text(encoding="utf-8") or "{}")
    else:
        data = {}

    entry = dict(station_diag)
    if extra:
        entry.update(extra)

    data.setdefault(param_name, [])
    data[param_name].append(entry)
    log_json_path.write_text(json.dumps(data, indent=2), encoding="utf-8")

    with open(log_txt_path, "a", encoding="utf-8") as f:
        f.write(f"[{param_name}] x0={entry.get('x0'):.6g}, dx={entry.get('dx_used'):.6g}, n={entry.get('n_points')}, "
                f"mean|x-x0|={entry.get('mean_abs_dx'):.6g}, std={entry.get('std_abs_dx'):.6g}, "
                f"r_min={entry.get('r_min'):.6g}, r_max={entry.get('r_max'):.6g}\n")

def title_from_template(template: str, var: str, param: str, x_in: float) -> str:
    try:
        return template.format(var=var, param=param, x_in=x_in)
    except Exception:
        return f"{var}: {param} — x={x_in:.2f} in"

def infer_param_name(path_like: Path) -> str:
    p = Path(path_like)
    return p.name
