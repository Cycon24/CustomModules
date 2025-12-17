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
RESIDUALS_HEADERS = ['Time_Iter', 'Outer_Iter', 'Inner_Iter', 'rms[Rho]',
       'rms[RhoU]', 'rms[RhoV]','rms[RhoW]', 'rms[RhoE]', 'rms[nu]']

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
    

    # # Read raw text & normalize
    # try:
    #     raw = src.read_text(encoding="utf-8-sig", errors="replace")
    # except Exception:
    #     raw = src.read_text(encoding="latin-1", errors="replace")
    # raw = raw.replace("\r\n", "\n").replace("\r", "\n")

    # # Strip trailing commas (phantom empty fields)
    # cleaned_lines = [re.sub(r",\s*$", "", line) for line in raw.split("\n")]
    # cleaned = "\n".join(cleaned_lines)

    # # Try 1: autodetect separator
    # try:
    #     df = pd.read_csv(
    #         io.StringIO(cleaned),
    #         sep=None, engine="python",
    #         skip_blank_lines=True,
    #         on_bad_lines="skip",
    #     )
    # except Exception as e1:
    #     print(f"[INFO] Auto-sep parse failed ({e1}); trying comma-or-whitespace regex...")
    #     # Try 2: comma OR whitespace
    #     try:
    #         df = pd.read_csv(
    #             io.StringIO(cleaned),
    #             sep=r",\s*|\s+",
    #             engine="python",
    #             skip_blank_lines=True,
    #             on_bad_lines="skip",
    #         )
    #     except Exception as e2:
    #         print(f"[WARN] Failed reading {src}: {e2}")
    #         return None

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

    # iter_col = _detect_iteration_column(df, user_hint=iter_col_hint)
    # if iter_col is None:
    #     print(f"[WARN] Could not detect iteration column in {param_dir / HISTORY_FILENAME}. Skipping residuals plot.")
    #     return

    # rms_cols = _find_rms_columns(df, iter_col)
    # if not rms_cols:
    #     print(f"[WARN] No residual columns found in {param_dir / HISTORY_FILENAME}. Skipping residuals plot.")
    #     return

    # use = df[[iter_col] + rms_cols].copy()
    # use[iter_col] = pd.to_numeric(use[iter_col], errors="coerce")
    # for c in rms_cols:
    #     use[c] = pd.to_numeric(use[c], errors="coerce")
    # use = use.dropna().sort_values(iter_col, kind="mergesort")
    # if use.empty:
    #     print(f"[WARN] No valid rows in {param_dir / HISTORY_FILENAME} after cleaning.")
    #     return

    plots_dir = _ensure_plots_dir(param_dir)
    fig, ax = plt.subplots()
    
    
    for c in df.columns.tolist():
        if c.startswith('rms'):
            # series = use[c]
            # yvals = np.power(10.0, series.to_numpy(float)) if _is_log10_series(series) else series.to_numpy(float)
            # xvals = use[iter_col].to_numpy(float)
            ax.plot(df["Inner_Iter"], df[c], label=_nice_residual_label(c))

    
    ax.set_title("Residuals")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("RMS")
    # ax.set_yscale("log")  # standard for residuals
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.set_ylim([-10,2])
    ax.minorticks_on()
    ax.legend()

    fig.tight_layout()
    outpath = plots_dir / f"Residuals.{img_format.lower()}"
    fig.savefig(outpath, dpi=dpi, format=img_format.lower())
    plt.close(fig)
    print(f"[OK] {param_dir.name}: saved {outpath.name}")
    
if __name__=="__main__":
    plot_residuals_history(Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\RANS_PERIODIC"), "png", 250)