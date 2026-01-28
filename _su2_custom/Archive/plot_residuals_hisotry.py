# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 12:30:17 2025

@author: BriceM
"""
from pathlib import Path
from typing import Optional
import re 


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HISTORY_FILENAME = "history.csv"
# RESIDUALS_HEADERS = ['Time_Iter', 'Outer_Iter', 'Inner_Iter', 'rms[Rho]',
#        'rms[RhoU]', 'rms[RhoV]','rms[RhoW]', 'rms[RhoE]', 'rms[nu]']
RESIDUALS_HEADERS = ['Inner_Iter', "Time(sec)", 'rms[Rho]',
       'rms[RhoU]', 'rms[RhoV]','rms[RhoW]', 'rms[RhoE]', 'rms[nu]', "CD", "CL"]

# def _ensure_plots_dir(param_dir: Path) -> Path:
#     plots_dir = param_dir / "Plots"
#     plots_dir.mkdir(exist_ok=True)
#     return plots_dir





def _find_rms_columns(df: pd.DataFrame, iter_col: str) -> list[str]:
    # residual cols often look like rms[...]
    cols = []
    for c in df.columns:
        if c == iter_col:
            continue
        cl = c.lower()
        if "rms[" in cl or cl.startswith("rms"):
            cols.append(c)
    return cols

def _is_log10_series(s: pd.Series) -> bool:
    s = pd.to_numeric(s, errors="coerce").dropna()
    return (not s.empty) and ((s < 0).mean() >= 0.2)




    
if __name__=="__main__":
    # plot_residuals_history(Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\BackPressure_SI_02_mdot"), "png", 300)
    plot_aero_coefficients(Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\BackPressure_SI_02_mdot"), "png", 300)