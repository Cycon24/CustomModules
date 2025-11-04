# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 11:23:45 2025

@author: BriceM
"""
from pathlib import Path
import shutil

def copy_mesh(filelocation, filename, newfilelocation, newfilename):
    """
    Copy a mesh pair (.msh and .su2) from one location to another with a new base name.

    Parameters
    ----------
    filelocation : str or Path
        Directory where the original mesh files live.
    filename : str
        Base filename (without extension) of the original mesh files.
        Expects: <filename>.msh and <filename>.su2 in filelocation.
    newfilelocation : str or Path
        Directory where the copied mesh files should be placed.
    newfilename : str
        New base filename (without extension) for the copied mesh files.
        Will create: <newfilename>.msh and <newfilename>.su2 in newfilelocation.

    Raises
    ------
    FileNotFoundError
        If either the .msh or .su2 source file does not exist.
    """
    src_dir = Path(filelocation)
    dst_dir = Path(newfilelocation)

    # Ensure destination directory exists
    dst_dir.mkdir(parents=True, exist_ok=True)

    exts = [".msh", ".su2"]
    for ext in exts:
        src = src_dir / f"{filename}{ext}"
        if not src.exists():
            raise FileNotFoundError(f"Source mesh file not found: {src}")

        dst = dst_dir / f"{newfilename}{ext}"
        shutil.copy2(src, dst)
    