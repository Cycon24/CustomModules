import argparse, yaml
from pathlib import Path
import pandas as pd

from cfd3d_utils import (
    read_surface_csv, add_cylindrical_about_x, extract_slab,
    ensure_dir, append_summary, infer_param_name
)

def main():
    ap = argparse.ArgumentParser(description="Extract slabbed cross-sections (probe surfaces) from entire_surface_restart.csv")
    ap.add_argument("--config", required=True, help="Path to post3d_config.yaml")
    ap.add_argument("--param_dir", default=".", help="Parameter folder containing entire_surface_restart.csv")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    param_dir = Path(args.param_dir)
    param_name = infer_param_name(param_dir)

    csv_path = param_dir / cfg["input"]["surface_csv"]
    columns = cfg["input"]["columns"]

    df = read_surface_csv(csv_path, columns)
    df = add_cylindrical_about_x(df)

    r_hub = float(cfg["geometry"]["r_hub"])
    r_pipe = float(cfg["geometry"]["r_pipe"])
    x_stations = list(cfg["slices"]["x_stations_in"])
    dx_tol = cfg["slices"].get("dx_tolerance_in", None)
    min_pts = int(cfg["slices"].get("min_points_per_slice", 3000))
    save_slice_csv = bool(cfg["slices"].get("save_slice_csv", True))

    probes_dir = param_dir / "Probes"
    logs_dir = param_dir / "Logs"
    ensure_dir(probes_dir)
    ensure_dir(logs_dir)

    for x0 in x_stations:
        slice_df, diags = extract_slab(df, float(x0), None if dx_tol is None else float(dx_tol), min_pts)
        if save_slice_csv:
            fn = probes_dir / f"slice_x={x0:.6f}_dx={diags['dx_used']:.6f}.csv"
            slice_df.to_csv(fn, index=False)

        append_summary(
            log_json_path=logs_dir / "post_summary.json",
            log_txt_path=logs_dir / "post_summary.txt",
            param_name=param_name,
            station_diag=diags,
            extra={
                "r_hub": r_hub, "r_pipe": r_pipe,
                "coverage_in_bounds": float(((slice_df['r']>=r_hub)&(slice_df['r']<=r_pipe)).mean()) if len(slice_df)>0 else 0.0,
            }
        )

if __name__ == "__main__":
    main()
