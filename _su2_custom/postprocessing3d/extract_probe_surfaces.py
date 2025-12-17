from pathlib import Path
from typing import Optional

import pandas as pd

from cfd3d_utils import (
    read_surface_csv,
    add_cylindrical_about_x,
    extract_slab,
    ensure_dir,
    append_summary,
    infer_param_name,
)


def extract_probe_surfaces(
    param_dir: Path,
    cfg: dict,
) -> None:
    """
    Generate slabbed cross-section CSVs ("probe surfaces") for each requested
    axial station in this param directory.

    param_dir : Path to the sweep param folder (contains entire_surface_restart.csv)
    cfg       : dict loaded from post3d_config.yaml
    """
    param_dir = Path(param_dir)
    param_name = infer_param_name(param_dir)

    # --- config unpack ---
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
        slice_df, diags = extract_slab(
            df=df,
            x0=float(x0),
            dx_tolerance=None if dx_tol is None else float(dx_tol),
            min_points=min_pts,
        )

        # save the slice data for reuse
        if save_slice_csv:
            fn = probes_dir / f"slice_x={x0:.6f}.csv"
            slice_df.to_csv(fn, index=False)

        # log diagnostics for summary
        coverage = (
            float(((slice_df["r"] >= r_hub) & (slice_df["r"] <= r_pipe)).mean())
            if len(slice_df) > 0
            else 0.0
        )

        append_summary(
            log_json_path=logs_dir / "post_summary.json",
            log_txt_path=logs_dir / "post_summary.txt",
            param_name=param_name,
            station_diag=diags,
            extra={
                "r_hub": r_hub,
                "r_pipe": r_pipe,
                "coverage_in_bounds": coverage,
            },
        )


