from pathlib import Path
from typing import Optional

from cfd3d_utils import load_main_config, load_units_config
from extract_probe_surfaces import extract_probe_surfaces
from plot_cross_section_contours import plot_cross_section_contours
from plot_radial_averages import plot_radial_averages


def run_full_postprocess(
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
    extract_probe_surfaces(param_dir, cfg)
    

    # stage 2: contours
    plot_cross_section_contours(param_dir, cfg, units_cfg)

    # stage 3: radial averages
    plot_radial_averages(param_dir, cfg, units_cfg)


if __name__=="__main__":
    run_full_postprocess(param_dir="C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\3D_Tests\\Pb_sweep\\Pb_1693",
                         cfg_path= r'C:\Users\BriceM\Documents\Modules\_su2_custom\3dpostprocessing\post3d_config.yaml')