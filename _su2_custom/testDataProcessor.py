# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 11:50:38 2025

@author: BriceM
"""
import postprocessing3d.Data3D_Processor as DP3
import postprocessing3d.Data3D_Plotter_Params as DPP3
import postprocessing3d.Data3D_Plotter_Sweep as DPS3
from pathlib import Path
cfg_path       = Path('C:\\Users\\BriceM\\Documents\\Modules\\_su2_custom\\3dpostprocessing\\post3d_config.yaml')
units_path     = r'C:\Users\BriceM\Documents\Modules\_su2_custom\3dpostprocessing\plot_units.yaml'
dpi = 250
sweep_path = r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\AoA_rt_sweep"
DP3.main(sweep_path, cfg_path, units_path=units_path)
DPP3.main(sweep_path, cfg_path, units_path=units_path, img_format="png", dpi=dpi)
DPS3.main(sweep_path, cfg_path, units_path=units_path, img_format="png", dpi=dpi)