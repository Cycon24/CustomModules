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

def _ensure_plots_dir(param_dir: Path) -> Path:
    plots_dir = param_dir / "Plots"
    plots_dir.mkdir(exist_ok=True)
    return plots_dir

def _sanitize_colname(name: str) -> str:
    s = str(name).strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in {"'", '"'}:
        s = s[1:-1]
    return s.strip()

def _read_history_csv(param_dir: Path) -> pd.DataFrame | None:
    """
    Read and sanitize history.csv robustly:
    - Normalize newlines and strip trailing commas that cause extra fields
    - Try auto-sep first; fallback to combined regex (comma OR whitespace)
    - Skip malformed lines rather than failing
    - Sanitize headers (strip spaces/quotes)
    """
    src = param_dir / HISTORY_FILENAME
    if not src.exists():
        print(f"[WARN] Missing {HISTORY_FILENAME} in {param_dir}. Skipping residuals plot.")
        return None
    
    df = pd.read_csv(src)
    df.columns = RESIDUALS_HEADERS[:len(df.columns.tolist())]
    

    if df is None or df.empty:
        print(f"[WARN] Empty/invalid {HISTORY_FILENAME} in {param_dir}.")
        return None

    # Clean headers
    # df = df.rename(columns=_sanitize_colname)
    return df

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

def _nice_residual_label(colname: str) -> str:
    m = re.search(r"rms\[(.+?)\]", colname, flags=re.IGNORECASE)
    return m.group(1) if m else colname

def plot_residuals_history(param_dir: Path, img_format: str, dpi: int, iter_col_hint: str | None = "Inner_Iter"):
    """
    Plot Residuals: RMS vs Iteration from history.csv
    - Iteration on x-axis (auto-picked or user hint)
    - Linear RMS on log y-axis
    """
    df = _read_history_csv(param_dir)
    if df is None:
        return

    plots_dir = _ensure_plots_dir(param_dir)
    fig, ax = plt.subplots()
    
    
    for c in df.columns.tolist():
        if c.startswith('rms'):
            ax.plot(df["Inner_Iter"], df[c], label=_nice_residual_label(c))

    ax.set_title("Residuals")
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
    
def plot_aero_coefficients(param_dir: Path, img_format: str, dpi: int):
    df = _read_history_csv(param_dir)
    if df is None:
        return
    
    

    plots_dir = _ensure_plots_dir(param_dir)
    fig, ax = plt.subplots()
    
    line_cl, = ax.plot(df["Inner_Iter"], df["CL"], label="CL",color='C0')
    ax.set_ylabel("Lift Coefficient - CL")
    
    ax2 = ax.twinx()
    line_cd, = ax2.plot(df["Inner_Iter"], df["CD"], label="CD", color="C1")
    ax2.set_ylabel("Drag Coefficient - CD")
    
    
    ax.set_title("Lift and Drag Coefficients")
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
    
if __name__=="__main__":
    # plot_residuals_history(Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\BackPressure_SI_02_mdot"), "png", 300)
    plot_aero_coefficients(Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\BackPressure_SI_02_mdot"), "png", 300)