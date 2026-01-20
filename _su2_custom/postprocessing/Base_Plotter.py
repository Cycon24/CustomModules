# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 10:08:39 2026

@author: BriceM
"""
from pathlib import Path
import matplotlib.pyplot as plt 
import pandas as pd


def plot_residuals(param_dir:Path, plt_res_params:dict):
    
    src = param_dir / plt_res_params["FileName"]
    
    if not src.exists():
        print(f"[WARN] Could not find residuals file in {param_dir}.")
        return None
    
    df = pd.read_csv(src)
    # df.columns = RESIDUALS_HEADERS[:len(df.columns.tolist())]