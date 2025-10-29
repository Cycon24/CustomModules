import argparse, yaml
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from cfd3d_utils import (
    ensure_dir, title_from_template, infer_param_name
)

def main():
    ap = argparse.ArgumentParser(description="Make contour plots on yz-slices for selected variables")
    ap.add_argument("--config", required=True, help="Path to post3d_config.yaml")
    ap.add_argument("--units", required=False, help="Path to plot_units.yaml (overrides config.plotting.units_yaml)")
    ap.add_argument("--param_dir", default=".", help="Parameter folder containing Probes/*.csv")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    units_yaml_path = Path(args.units) if args.units else Path(cfg["plotting"]["units_yaml"])
    units_cfg = yaml.safe_load(units_yaml_path.read_text(encoding="utf-8")) if units_yaml_path.exists() else {"title":{"template":"{var}: {param} — x={x_in:.2f} in"}, "units":{}, "limits":{}}

    param_dir = Path(args.param_dir)
    param_name = infer_param_name(param_dir)

    contours_dir = param_dir / (units_cfg.get("naming",{}).get("contours_dir","Plots/Contours"))
    ensure_dir(contours_dir)

    dpi = int(cfg["plotting"].get("dpi", 200))
    cmap = cfg["plotting"].get("cmap", "viridis")
    levels = cfg["contours"].get("levels", 50)
    variables = list(cfg["contours"]["variables"])

    probe_dir = param_dir / (units_cfg.get("naming",{}).get("probes_dir","Probes"))
    probe_files = sorted(probe_dir.glob("slice_x=*_dx=*.csv"))
    if not probe_files:
        raise FileNotFoundError(f"No probe CSVs found in {probe_dir}")

    title_tmpl = units_cfg.get("title",{}).get("template", "{var}: {param} — x={x_in:.2f} in")
    yz_limits = units_cfg.get("limits",{}).get("yz", {})
    y_lim = yz_limits.get("y_in", None)
    z_lim = yz_limits.get("z_in", None)

    for pf in probe_files:
        name = pf.stem
        try:
            x_in = float(name.split("slice_x=")[1].split("_dx=")[0])
        except Exception:
            x_in = np.nan

        df = pd.read_csv(pf)
        tri = mtri.Triangulation(df["y"].to_numpy(), df["z"].to_numpy())

        for var in variables:
            if var not in df.columns:
                alias = {
                    "Vmag":"Vmag","Vx":"Vx","Vr":"Vr","Vtheta":"Vtheta",
                    "Mach":"Mach","Pressure":"Pressure","Temperature":"Temperature",
                    "Momentum_x":"Momentum_x","Momentum_r":"Momentum_r","Momentum_theta":"Momentum_theta",
                    "alpha_deg":"alpha_deg"
                }.get(var, None)
                if alias is None or alias not in df.columns:
                    continue

            vals = df[var].to_numpy()
            fig, ax = plt.subplots(figsize=(6,5), dpi=dpi)
            cs = ax.tricontourf(tri, vals, levels=levels, cmap=cmap)

            vlims = units_cfg.get("limits",{}).get(var, None)
            if vlims and isinstance(vlims, (list, tuple)) and len(vlims)==2:
                cs.set_clim(vlims[0], vlims[1])

            cbar = fig.colorbar(cs, ax=ax)
            unit = units_cfg.get("units",{}).get(var, "")
            if unit:
                cbar.set_label(unit)

            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("y [in]")
            ax.set_ylabel("z [in]")
            if y_lim: ax.set_xlim(y_lim)
            if z_lim: ax.set_ylim(z_lim)

            title = title_from_template(title_tmpl, var=var, param=param_name, x_in=x_in)
            ax.set_title(title)

            out_dir = contours_dir / f"x={x_in:.2f}"
            ensure_dir(out_dir)
            out_path = out_dir / f"{var}.png"
            fig.savefig(out_path, bbox_inches="tight")
            plt.close(fig)

if __name__ == "__main__":
    main()
