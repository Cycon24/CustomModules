from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from cfd3d_utils import (
    ensure_dir,
    radial_bin_stats,
    infer_param_name,
    load_units_config,
)


def plot_radial_averages(
    param_dir: Path,
    cfg: dict,
    units_cfg: Optional[dict] = None,
) -> None:
    """
    Build radially-averaged profiles for each saved probe surface.
    One figure per variable, with multiple station curves on the same axes.
    Also writes CSV stacks of the binned stats.

    param_dir : sweep param folder
    cfg       : main config dict
    units_cfg : plot_units.yaml dict, or None to auto-load
    """
    param_dir = Path(param_dir)
    param_name = infer_param_name(param_dir)

    # resolve units_cfg if not passed in
    if units_cfg is None:
        units_cfg = load_units_config(
            units_path=None,
            fallback_from_cfg=Path(cfg["plotting"]["units_yaml"]) if "plotting" in cfg else None,
        )

    radial_dir = param_dir / units_cfg.get("naming", {}).get("radial_dir", "Plots/Radial")
    ensure_dir(radial_dir)

    dpi = int(cfg["plotting"].get("dpi", 200))
    variables = list(cfg["radial_averages"]["variables"])
    Nr = int(cfg["radial_averages"]["Nr_bins"])

    r_hub = float(cfg["geometry"]["r_hub"])
    r_pipe = float(cfg["geometry"]["r_pipe"])

    probe_dir = param_dir / units_cfg.get("naming", {}).get("probes_dir", "Probes")
    probe_files = sorted(probe_dir.glob("slice_x=*.csv"))
    if not probe_files:
        raise FileNotFoundError(f"No probe CSVs found in {probe_dir}")

    # build radial stats for each station
    per_station = {}
    stations = []
    for pf in probe_files:
        name = pf.stem  # slice_x=..._dx=...
        try:
            x_in = float(name.split("slice_x=")[1].split("_dx=")[0])
        except Exception:
            x_in = np.nan

        df = pd.read_csv(pf)
        stats = radial_bin_stats(df, r_hub, r_pipe, Nr, variables)
        per_station[x_in] = stats
        stations.append(x_in)

    stations = sorted([s for s in stations if not np.isnan(s)])
    if not stations:
        raise RuntimeError("No valid x-station values parsed from probe file names.")

    # reference radius axis for plotting
    r_mid = per_station[stations[0]]["r_mid"].to_numpy()

    for var in variables:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)

        for x_in in stations:
            stats = per_station[x_in]
            yvals = stats[f"{var}_mean"].to_numpy()
            ax.plot(r_mid, yvals, label=f"x={x_in:.2f} in")

        unit_lbl = units_cfg.get("units", {}).get(var, "")
        ax.set_xlabel("r [in]")
        ax.set_ylabel(f"{var} {'[' + unit_lbl + ']' if unit_lbl else ''}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

        # multi-station title template
        tmpl = units_cfg.get("title", {}).get("template_multi", "{var}: {param}")
        title = tmpl.format(var=var, param=param_name)
        ax.set_title(title)

        out_path = radial_dir / f"{var}.png"
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)

        # also save combined CSV across all stations for this var
        rows = []
        for x_in in stations:
            stats = per_station[x_in]
            dfv = stats[["r_mid", f"{var}_mean", f"{var}_std", "count"]].copy()
            dfv["x_in"] = x_in
            rows.append(dfv)
        stack = pd.concat(rows, ignore_index=True)
        stack.to_csv(radial_dir / f"{var}_radial_stats.csv", index=False)
