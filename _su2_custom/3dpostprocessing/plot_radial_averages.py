import argparse, yaml
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from cfd3d_utils import (
    ensure_dir, radial_bin_stats, infer_param_name
)

def main():
    ap = argparse.ArgumentParser(description="Radially averaged line plots per plane; one figure per variable with one curve per station")
    ap.add_argument("--config", required=True, help="Path to post3d_config.yaml")
    ap.add_argument("--units", required=False, help="Path to plot_units.yaml (overrides config.plotting.units_yaml)")
    ap.add_argument("--param_dir", default=".", help="Parameter folder containing Probes/*.csv")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    units_yaml_path = Path(args.units) if args.units else Path(cfg["plotting"]["units_yaml"])
    units_cfg = yaml.safe_load(units_yaml_path.read_text(encoding="utf-8")) if units_yaml_path.exists() else {"title":{"template_multi":"{var}: {param}"}}

    param_dir = Path(args.param_dir)
    param_name = infer_param_name(param_dir)

    radial_dir = param_dir / (units_cfg.get("naming",{}).get("radial_dir","Plots/Radial"))
    ensure_dir(radial_dir)

    dpi = int(cfg["plotting"].get("dpi", 200))
    variables = list(cfg["radial_averages"]["variables"])
    Nr = int(cfg["radial_averages"]["Nr_bins"])
    r_hub = float(cfg["geometry"]["r_hub"])
    r_pipe = float(cfg["geometry"]["r_pipe"])

    probe_dir = param_dir / (units_cfg.get("naming",{}).get("probes_dir","Probes"))
    probe_files = sorted(probe_dir.glob("slice_x=*_dx=*.csv"))
    if not probe_files:
        raise FileNotFoundError(f"No probe CSVs found in {probe_dir}")

    per_station = {}
    stations = []
    for pf in probe_files:
        name = pf.stem
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
        raise RuntimeError("No valid stations parsed from probe file names.")

    r_mid = per_station[stations[0]]["r_mid"].to_numpy()

    for var in variables:
        fig, ax = plt.subplots(figsize=(6,4), dpi=dpi)
        for x_in in stations:
            stats = per_station[x_in]
            y = stats[f"{var}_mean"].to_numpy()
            ax.plot(r_mid, y, label=f"x={x_in:.2f} in")
        unit = units_cfg.get("units",{}).get(var, "")
        ax.set_xlabel("r [in]")
        ax.set_ylabel(f"{var} {f'[{unit}]' if unit else ''}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)

        tmpl = units_cfg.get("title",{}).get("template_multi", "{var}: {param}")
        title = tmpl.format(var=var, param=param_name)
        ax.set_title(title)

        out_path = radial_dir / f"{var}.png"
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)

        rows = []
        for x_in in stations:
            stats = per_station[x_in]
            dfv = stats[["r_mid", f"{var}_mean", f"{var}_std", "count"]].copy()
            dfv["x_in"] = x_in
            rows.append(dfv)
        stack = pd.concat(rows, ignore_index=True)
        stack.to_csv(radial_dir / f"{var}_radial_stats.csv", index=False)

if __name__ == "__main__":
    main()
