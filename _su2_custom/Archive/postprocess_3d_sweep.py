from pathlib import Path
from typing import Optional

if __name__=="__main__":
    from cfd3d_utils import load_main_config, load_units_config
    from extract_probe_surfaces import extract_probe_surfaces
    from plot_cross_section_contours import plot_cross_section_contours
    from plot_radial_averages import plot_radial_averages
    from plot_residuals_hisotry import plot_residuals_history
else:
    from postprocessing3d.cfd3d_utils import load_main_config, load_units_config
    from postprocessing3d.extract_probe_surfaces import extract_probe_surfaces
    from postprocessing3d.plot_cross_section_contours import plot_cross_section_contours
    from postprocessing3d.plot_radial_averages import plot_radial_averages


def run_postprocess_param(
    param_dir: Path,
    cfg_path: Path,
    units_path: Optional[Path] = None,
) -> None:
    """
    One-call convenience function to:
      1. slice & save probe surfaces
      2. generate contour plots at each station
      3. generate radial-average plots

    param_dir  : path to e.g. "AoA_10" directory containing entire_surface_restart.csv
    cfg_path   : path to post3d_config.yaml
    units_path : optional override path to plot_units.yaml
    """
    param_dir = Path(param_dir)
    cfg_path = Path(cfg_path)

    # load configs
    cfg = load_main_config(cfg_path)
    units_cfg = load_units_config(
        units_path=units_path,
        fallback_from_cfg=Path(cfg["plotting"]["units_yaml"]) if "plotting" in cfg else None,
    )

    # stage 1: generate probe surfaces & summary logs
    # extract_probe_surfaces(param_dir, cfg)
    

    # # stage 2: contours
    # plot_cross_section_contours(param_dir, cfg, units_cfg)

    # # stage 3: radial averages
    # plot_radial_averages(param_dir, cfg, units_cfg)
    
    # stage 4: residuals
    plot_residuals_history(param_dir, "png", dpi=250)
    
def _iter_param_dirs(sweep_root: Path):
    """Yield immediate subdirectories (Param folders)."""
    for child in sorted(sweep_root.iterdir()):
        if child.is_dir():
            yield child
    
def run_postprocess_sweep(sweep_dir: Path,
    cfg_path: Path,
    units_path: Optional[Path] = None,
) -> None:
    
    sr_path = Path(sweep_dir).resolve()
    if not sr_path.exists():
        raise FileNotFoundError(f"[ERROR] sweep_root not found: {sr_path}")

    any_found = False
    for param_dir in _iter_param_dirs(sr_path):
        any_found = True
        try:
            run_postprocess_param(
                param_dir=param_dir,
                cfg_path=cfg_path,
                units_path=units_path)
        except FileNotFoundError as e:
            raise Warning(f"No file or directory: {e}")

if __name__=="__main__":
    run_postprocess_sweep(sweep_dir="C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\3D_Tests\\AoA_rt_sweep",
                         cfg_path= 'post3d_config.yaml')