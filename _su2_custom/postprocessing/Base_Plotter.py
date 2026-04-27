# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 10:08:39 2026

@author: BriceM
"""
from pathlib import Path
import matplotlib.pyplot as plt 
import matplotlib.tri as mtri
import numpy as np
import pandas as pd
import re
from typing import Optional, Dict, Any

# Local
import processing_tools as proc 
import process_mass_flow as proc_flow

CFG_DEF_NAME = "post3d_config.yaml"

# =============================================================================
# Processing
# =============================================================================
def extract_probe_surfaces(
    param_dir: Path,
    cfg: dict | None,
    cfg_path: Path | None = None
    ) -> None:
    """
    Generate slabbed cross-section CSVs ("probe surfaces") for each requested
    axial station in this param directory.

    param_dir : Path to the sweep param folder (contains entire_surface_restart.csv)
    cfg       : dict loaded from post3d_config.yaml
    """
    # Load cfg file if dict isnt provided
    if cfg == None:
        if cfg_path == None:
            # Get the absolute path of the current file
            current_file_path = Path(__file__).resolve()
            # If you only need the directory containing the file
            cfg_path = current_file_path.parent.parent / "config" / CFG_DEF_NAME
        
        cfg = proc.load_main_config(cfg_path)
    
    
    param_dir = Path(param_dir)
    param_name = proc.infer_param_name(param_dir)

    # --- config unpack ---
    csv_path = param_dir / cfg["input"]["surface_csv"]
    columns = cfg["input"]["columns"]
    scale_factor = cfg["units"]["scale_factor"]
    units = cfg["units"].get("final_units", "")
    
    df = proc.read_surface_csv(csv_path, columns, scale_factor)
    df = proc.add_cylindrical_about_x(df)

    r_hub = float(cfg["geometry"]["r_hub"])
    r_pipe = float(cfg["geometry"]["r_pipe"])

    x_stations_in = list(cfg["slices"]["x_stations_in"])
    dx_tol = cfg["slices"].get("dx_tolerance", None)
    min_pts = int(cfg["slices"].get("min_points_per_slice", 3000))
    save_slice_csv = bool(cfg["slices"].get("save_slice_csv", True))
    

    probes_dir = param_dir / "Probes"
    logs_dir = param_dir / "Logs"
    proc.ensure_dir(probes_dir)
    proc.ensure_dir(logs_dir)

    for x0 in x_stations_in:
        slice_df, diags = proc.extract_slab(
            df=df,
            x0=float(x0),
            dx_tolerance=None if dx_tol is None else float(dx_tol),
            min_points=min_pts,
        )

        # save the slice data for reuse
        if save_slice_csv:
            fn = probes_dir / f"slice_x={x0:.6f}_{units}.csv"
            slice_df.to_csv(fn, index=False)

        # log diagnostics for summary
        coverage = (
            float(((slice_df["r"] >= r_hub) & (slice_df["r"] <= r_pipe)).mean())
            if len(slice_df) > 0
            else 0.0
        )

        proc.append_summary(
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
    return None
# =============================================================================
# Plotting
# =============================================================================
def plot_residuals_history(param_dir: Path, 
                           img_format: str = "png", 
                           dpi: int = 300, 
                           iter_col_hint: str | None = "Inner_Iter"):
    """
    Plot Residuals: RMS vs Iteration from history.csv
    - Iteration on x-axis (auto-picked or user hint)
    - Linear RMS on log y-axis
    """
    df = proc.read_history_csv(param_dir)
    if df is None:
        print(f"[Warn] {param_dir.name} could not be imported.")
        return

    plots_dir = proc.ensure_dir(param_dir / "Plots")
    fig, ax = plt.subplots()
    
    
    def _nice_residual_label(colname: str) -> str:
        m = re.search(r"rms\[(.+?)\]", colname, flags=re.IGNORECASE)
        return m.group(1) if m else colname

    for c in df.columns.tolist():
        if c.startswith('rms'):
            ax.plot(df["Inner_Iter"], df[c], label=_nice_residual_label(c))

    ax.set_title(f"Residuals\n {param_dir.name}")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("log(RMS)")
    plt.minorticks_on()
    plt.grid(visible=True, which='major',color='k', linestyle='-')
    plt.grid(visible=True, which='minor', linestyle='--', linewidth=0.5)
    ax.set_ylim([-10,2])
    ax.legend(ncol=2)

    fig.tight_layout()
    outpath = plots_dir / f"Residuals.{img_format.lower()}"
    fig.savefig(outpath, dpi=dpi, format=img_format.lower())
    plt.close(fig)
    print(f"[OK] {param_dir.name}: saved {outpath.name}")
    
def plot_aero_coefficients(param_dir: Path, 
                           img_format: str= "png", 
                           dpi: int = 300):
    
    df = proc.read_history_csv(param_dir)
    if df is None:
        print(f"[Warn] {param_dir.name} could not be imported.")
        return
    

    plots_dir = proc.ensure_dir(param_dir / "Plots")
    fig, ax = plt.subplots()
    
    line_cl, = ax.plot(df["Inner_Iter"], df["CL"], label="CL",color='C0')
    ax.set_ylabel("Lift Coefficient - CL")
    
    ax2 = ax.twinx()
    line_cd, = ax2.plot(df["Inner_Iter"], df["CD"], label="CD", color="C1")
    ax2.set_ylabel("Drag Coefficient - CD")
    
    
    ax.set_title(f"Lift and Drag Coefficients\n {param_dir.name}")
    ax.set_xlabel("Iteration")
    
    ax.minorticks_on()
    ax.grid(visible=True, which='major',color='k', linestyle='-')
    ax.grid(visible=True, which='minor', linestyle='--', linewidth=0.5)
    # ax.set_ylim([-10,2])
    # Combined legend
    lines = [line_cl, line_cd]
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, ncol=2)

    fig.tight_layout()
    outpath = plots_dir / f"CoefficientResiduals.{img_format.lower()}"
    fig.savefig(outpath, dpi=dpi, format=img_format.lower())
    plt.close(fig)
    print(f"[OK] {param_dir.name}: saved {outpath.name}")
    
    
def plot_cross_section_contours(
    param_dir: Path,
    cfg: dict | None,
    units_cfg: Optional[dict] = None,
    cfg_path: Path | None = None  
) -> None:
    '''
    For each saved probe surface in param_dir/Probes, make (y,z) tricontourf plots
    for each variable in cfg["contours"]["variables"].

    Parameters
    ----------
    param_dir : Path
        Sweep param folder, or folder that contains the immediate CFD data.
    cfg : dict | None
        Plotting config dict.
    units_cfg : Optional[dict], optional
        Dict from plot_units.yaml (or None to auto-load using cfg["plotting"]["units_yaml"]). The default is None.
    cfg_path : Path | None, optional
        Path to the plot config file, must include the filename and extension "name.yaml". 
        The default is None which will utilize the internal package file location.

    Raises
    ------
    FileNotFoundError
        Raised when a CSV file cannot be found.

    Returns
    -------
    None
       
    '''
   
    param_dir = Path(param_dir)
    param_name = proc.infer_param_name(param_dir)
    
    # Load cfg file if dict isnt provided
    if cfg == None:
        if cfg_path == None:
            # Get the absolute path of the current file
            current_file_path = Path(__file__).resolve()
            # If you only need the directory containing the file
            cfg_path = current_file_path.parent.parent / "config" / CFG_DEF_NAME
        
        cfg = proc.load_main_config(cfg_path)

    # resolve units_cfg if not passed
    if units_cfg is None:
        units_cfg = proc.load_units_config(
            units_path=None,
            fallback_from_cfg=Path(cfg["plotting"]["units_yaml"]) if "plotting" in cfg else None,
        )

    contours_dir = param_dir / units_cfg.get("naming", {}).get("contours_dir", "Plots/Contours")
    proc.ensure_dir(contours_dir)

    dpi = int(cfg["plotting"].get("dpi", 200))
    cmap = cfg["plotting"].get("cmap", "viridis")
    levels = cfg["contours"].get("levels", 50)
    variables = list(cfg["contours"]["variables"])

    probe_dir = param_dir / units_cfg.get("naming", {}).get("probes_dir", "Probes")
    probe_files = sorted(probe_dir.glob("slice_x=*.csv"))
    if not probe_files:
        raise FileNotFoundError(f"No probe CSVs found in {probe_dir}")

    title_tmpl = units_cfg.get("title", {}).get("template", "{var}: {param} \n x={x_in:.2f} in")
    yz_limits = units_cfg.get("limits", {}).get("yz", {})
    y_lim = yz_limits.get("y_in", None)
    z_lim = yz_limits.get("z_in", None)

    for pf in probe_files:
        name = pf.stem  # "slice_x=<x>_dx=<dx>"
       
        try:
            x_in = float(name.split("slice_x=")[1].split("_")[0])
        except Exception:
            x_in = np.nan
            
        print(f"[INFO] Contour Plotter: plotting {param_name} slice x={x_in}")
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

            title = proc.title_from_template(title_tmpl, var=var, param=param_name, x_in=x_in)
            ax.set_title(title)

            out_dir = contours_dir / f"x={x_in:.2f}"
            proc.ensure_dir(out_dir)

            out_path = out_dir / f"{var}.png"
            fig.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
            
            
def plot_radial_averages(
    param_dir: Path,
    cfg: dict | None,
    units_cfg: Optional[dict] = None,
    cfg_path: Path | None = None  
) -> None:
    '''
    Build radially-averaged profiles for each saved probe surface (stations).
    One figure per variable denoted in the config file, with eacg station curves on the same plot.
    Also writes CSV stacks of the binned stats.

    Parameters
    ----------
    param_dir : Path
        Sweep param folder, or folder that contains the immediate CFD data.
    cfg : dict | None
        Plotting config dict.
    units_cfg : Optional[dict], optional
        Dict from plot_units.yaml (or None to auto-load using cfg["plotting"]["units_yaml"]). The default is None.
    cfg_path : Path | None, optional
        Path to the plot config file, must include the filename and extension "name.yaml". 
        The default is None which will utilize the internal package file location.

    Raises
    ------
    FileNotFoundError
        Raised when a CSV file cannot be found.

    Returns
    -------
    None
       
    ''' 
    param_dir = Path(param_dir)
    param_name = proc.infer_param_name(param_dir)
    
    # Load cfg file if dict isnt provided
    if cfg == None:
        if cfg_path == None:
            # Get the absolute path of the current file
            current_file_path = Path(__file__).resolve()
            # If you only need the directory containing the file
            cfg_path = current_file_path.parent.parent / "config" / CFG_DEF_NAME
        
        cfg = proc.load_main_config(cfg_path)

    # resolve units_cfg if not passed in
    if units_cfg is None:
        units_cfg = proc.load_units_config(
            units_path=None,
            fallback_from_cfg=Path(cfg["plotting"]["units_yaml"]) if "plotting" in cfg else None,
        )

    radial_dir = param_dir / units_cfg.get("naming", {}).get("radial_dir", "Plots/Radial")
    proc.ensure_dir(radial_dir)

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
            x_in = float(name.split("slice_x=")[1].split("_")[0])
        except Exception:
            x_in = np.nan

        df = pd.read_csv(pf)
        stats = proc.radial_bin_stats(df, r_hub, r_pipe, Nr, variables)
        per_station[x_in] = stats
        stations.append(x_in)

    stations = sorted([s for s in stations if not np.isnan(s)])
    if not stations:
        raise RuntimeError("No valid x-station values parsed from probe file names.")

    # reference radius axis for plotting
    r_mid = per_station[stations[0]]["r_mid"].to_numpy()

    for var in variables:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)
        print(f"[INFO] Radial Plotter: plotting {param_name} var={var}")

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
    
# def main(
#     sweep_root: str,
#     cfg_path: str,
#     units_path: Optional[str] = None,
# ) -> None:
#     """
#     Perform the 3-D preprocessing stage across an entire sweep.

#     This is analogous to Data_Processor.main(...) in your 2-D pipeline.

#     For each param folder under sweep_root:
#       1. reads entire_surface_restart.csv
#       2. generates axial "probe surface" slabs at each requested x-station
#          (using dx_tolerance_in / min_points_per_slice from the config)
#       3. logs slice diagnostics in Logs/post_summary.*

#     Args
#     ----
#     sweep_root : str
#         Path to the sweep directory containing subfolders like AoA_10, AoA_12, ...
#     cfg_path : str
#         Path to post3d_config.yaml
#     units_path : Optional[str]
#         Currently unused in this stage, but kept for API symmetry with
#         Data3D_Plotter_Params.main and future expansion (for example, if we
#         someday want unit-aware filtering here).
#     """
#     sweep_root = Path(sweep_root).resolve()
#     if not sweep_root.exists():
#         raise FileNotFoundError(f"[ERROR] sweep_root not found: {sweep_root}")

#     cfg_path = Path(cfg_path).resolve()
#     if not cfg_path.exists():
#         raise FileNotFoundError(f"[ERROR] cfg_path not found: {cfg_path}")

#     cfg = load_main_config(cfg_path)

#     any_found = False
#     for param_dir in _iter_param_dirs(sweep_root):
#         any_found = True
#         _run_for_param(param_dir, cfg)

#     if not any_found:
#         warnings.warn(f"[WARN] No subfolders found under {sweep_root}")
#     else:
#         print("[INFO] 3D Processor complete.")

def process_flow_values(param_dir: Path,
                        print_info:bool = True,
                        save_info:bool = True):
    filename = "entire_surface_restart.csv"
    
    imported    = proc_flow.import_csv_to_df(filename, param_dir)
    inlet       = proc_flow.extract_points_in_plane(0, imported, tol=0.01e-3)
    outlet      = proc_flow.extract_points_in_plane(17*2.54/100, imported, tol=0.01e-3)
    TE_slice    = proc_flow.extract_points_in_plane(8 *2.54/100, imported, tol=0.10e-3)
    
    mdot, diag      = proc_flow.mass_flow_rate_yz(inlet, return_diagnostics=True, gc=1)
    mdot_o, diag_o  = proc_flow.mass_flow_rate_yz(outlet, return_diagnostics=True, gc=1)
    # R_refs = np.linspace(0.2*2.54/100, 2.0*2.54/100, 100)
    # SNs = []
    # for R_ref in R_refs:
    SN, diag_SN         = proc_flow.swirl_number_yz(outlet, return_diagnostics=True, R_ref=None) #2.0*2.54/100)
    SN_TE, diag_SN_TE   = proc_flow.swirl_number_yz(TE_slice, return_diagnostics=True, R_ref=None) #2.0*2.54/100)
   
        # SNs.append(SN)
    Ptots_i, inlet =  proc_flow.total_pressure_slice(inlet, 1.4)
    Ptots_o, outlet =  proc_flow.total_pressure_slice(outlet, 1.4)
    dPtots = Ptots_o - Ptots_i 
    pi_s = Ptots_o / Ptots_i
    pi_avg = np.average(pi_s)
    
    dmdot = 6*(mdot_o - mdot)
    # print(f"mdot_s = {mdot:.4f} lbm/s")
    # print(f"mdot_t = {6*mdot:.4f} lbm/s")
    
    
    lines = [
    r"- $\dot{m}$" + f" = {6*mdot:.4f} kg/s",
    r"- $\Delta\dot{m}$" + f" = {dmdot:.4f} kg/s",
    r"- $\pi_{avg}$" + f" = {pi_avg:.4f}",
    r"- $\pi_{min}$" + f" = {np.min(pi_s):.4f}",
    r"- $\pi_{max}$" + f" = {np.max(pi_s):.4f}",
    f"- SN     = {SN:.4f}",
    r"- $\alpha_{SN}$" + f"     = {proc_flow.swirl_angle(SN):.4f}" + r"$^o$",
    r"- $SN_{TE}$   = " + f"{SN_TE:.4f}",
    r"- $\alpha_{SN, TE}$" + f"     = {proc_flow.swirl_angle(SN_TE):.4f}" + r"$^o$"
    ]
    print(f"[Info] Flow Results calculated for {param_dir.name}")
    with open(param_dir / "Flow_Results.txt", "w") as f:
        for line in lines:
            if print_info: print(line)
            if save_info: f.write(line + "\n")
    
    return None

def process_single_run(param_dir: Path, 
                       cfg_path: Path | None = None,
                       units_path: Path | None = None,
                       ):
    if units_path != None:
        units_cfg = proc.load_units_config(units_path, fallback_from_cfg=True)
    else: 
        units_cfg = None
    
    print(f"[Info] Beginning processing for {param_dir.name}")
    process_flow_values(param_dir)
    extract_probe_surfaces(param_dir, cfg=None, cfg_path=cfg_path)
    plot_residuals_history(param_dir)
    plot_aero_coefficients(param_dir)
    plot_cross_section_contours(param_dir, cfg=None, units_cfg=units_cfg, cfg_path=cfg_path)
    plot_radial_averages(param_dir, cfg=None, units_cfg=units_cfg, cfg_path=cfg_path)
    return None




if __name__=="__main__":
        param_dir = Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\MinThickTests\Test06")
        # extract_probe_surfaces(param_dir, cfg=None)
        # plot_residuals_history(param_dir)
        # plot_aero_coefficients(param_dir)
        # plot_cross_section_contours(param_dir, cfg=None)
        # plot_radial_averages(param_dir, cfg=None)
        process_single_run(param_dir)
        