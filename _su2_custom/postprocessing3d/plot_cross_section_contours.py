from pathlib import Path
from typing import Optional, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from cfd3d_utils import (
    ensure_dir,
    title_from_template,
    infer_param_name,
    load_units_config,
)


def plot_cross_section_contours(
    param_dir: Path,
    cfg: dict,
    units_cfg: Optional[dict] = None,
) -> None:
    """
    For each saved probe surface in param_dir/Probes, make (y,z) tricontourf plots
    for each variable in cfg["contours"]["variables"].

    param_dir : sweep param folder
    cfg       : main config dict
    units_cfg : dict from plot_units.yaml (or None to auto-load using cfg["plotting"]["units_yaml"])
    """
    param_dir = Path(param_dir)
    param_name = infer_param_name(param_dir)

    # resolve units_cfg if not passed
    if units_cfg is None:
        units_cfg = load_units_config(
            units_path=None,
            fallback_from_cfg=Path(cfg["plotting"]["units_yaml"]) if "plotting" in cfg else None,
        )

    contours_dir = param_dir / units_cfg.get("naming", {}).get("contours_dir", "Plots/Contours")
    ensure_dir(contours_dir)

    dpi = int(cfg["plotting"].get("dpi", 200))
    cmap = cfg["plotting"].get("cmap", "viridis")
    levels = cfg["contours"].get("levels", 50)
    variables = list(cfg["contours"]["variables"])

    probe_dir = param_dir / units_cfg.get("naming", {}).get("probes_dir", "Probes")
    probe_files = sorted(probe_dir.glob("slice_x=*.csv"))
    if not probe_files:
        raise FileNotFoundError(f"No probe CSVs found in {probe_dir}")

    title_tmpl = units_cfg.get("title", {}).get("template", "{var}: {param} — x={x_in:.2f} in")
    yz_limits = units_cfg.get("limits", {}).get("yz", {})
    y_lim = yz_limits.get("y_in", None)
    z_lim = yz_limits.get("z_in", None)

    for pf in probe_files:
        name = pf.stem  # "slice_x=<x>_dx=<dx>"
        try:
            x_in = float(name.split("slice_x=")[1].split("_dx=")[0])
        except Exception:
            x_in = np.nan

        df = pd.read_csv(pf)
        tri = mtri.Triangulation(df["y"].to_numpy(), df["z"].to_numpy())

        for var in variables:
            # graceful skip if column missing
            if var not in df.columns:
                # allow aliases
                alias_map = {
                    "Vmag": "Vmag",
                    "Vx": "Vx",
                    "Vr": "Vr",
                    "Vtheta": "Vtheta",
                    "Mach": "Mach",
                    "Pressure": "Pressure",
                    "Temperature": "Temperature",
                    "Momentum_x": "Momentum_x",
                    "Momentum_r": "Momentum_r",
                    "Momentum_theta": "Momentum_theta",
                    "alpha_deg": "alpha_deg",
                }
                alias = alias_map.get(var, None)
                if alias is None or alias not in df.columns:
                    continue

            vals = df[var].to_numpy()
            fig, ax = plt.subplots(figsize=(6, 5), dpi=dpi)
            cs = ax.tricontourf(tri, vals, levels=levels, cmap=cmap)

            # optional color limits overrides
            vlims = units_cfg.get("limits", {}).get(var, None)
            if (
                vlims
                and isinstance(vlims, (list, tuple))
                and len(vlims) == 2
                and all(v is not None for v in vlims)
            ):
                cs.set_clim(vlims[0], vlims[1])

            cbar = fig.colorbar(cs, ax=ax)
            unit_lbl = units_cfg.get("units", {}).get(var, "")
            if unit_lbl:
                cbar.set_label(unit_lbl)

            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("y [in]")
            ax.set_ylabel("z [in]")

            if y_lim:
                ax.set_xlim(y_lim)
            if z_lim:
                ax.set_ylim(z_lim)

            title = title_from_template(title_tmpl, var=var, param=param_name, x_in=x_in)
            ax.set_title(title)

            out_dir = contours_dir / f"x={x_in:.2f}"
            ensure_dir(out_dir)

            out_path = out_dir / f"{var}.png"
            fig.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
