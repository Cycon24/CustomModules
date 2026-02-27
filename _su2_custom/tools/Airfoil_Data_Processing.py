# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 16:37:32 2026

@author: BriceM
"""
from pathlib import Path 
import pandas as pd
import matplotlib.pyplot as plt 
from typing import Dict, List, Tuple, Optional, Union
import numpy as np
'''
Maybe need to read in the .su2 and pull all points within the elements under marker=airfoils and utilize
those point IDs for the generation of a csv for the airfoil?
'''



def read_su2_markers_as_dfs(filepath: Union[str, Path]) -> List[pd.DataFrame]:
    """
    Parse an SU2 mesh file and return a list of DataFrames, one per MARKER_TAG.

    Each DataFrame:
      - df.name == MARKER_TAG string
      - columns: ["elem", "p1", "p2", "p3"]
        * elem: integer element type from SU2 marker element lines (e.g., 5 for TRIANGLE)
        * p1..p3: point indices defining the element
          - For 2-node elements (e.g., line): p3 is <NA>
          - For >3-node elements: only first 3 node indices are kept (others ignored)

    Notes:
      - This reads marker (boundary) elements under the MARKER_TAG blocks.
      - Assumes standard SU2 format: MARKER_TAG= <name>, MARKER_ELEMS= <n>, then n element lines.

    Parameters
    ----------
    filepath : str | Path
        Path to the .su2 mesh file.

    Returns
    -------
    List[pandas.DataFrame]
        List of marker DataFrames. The marker name is stored in df.name.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"SU2 mesh file not found: {filepath}")

    # Storage: marker_name -> list of rows (elem, p1, p2, p3)
    marker_rows: Dict[str, List[Tuple[int, Optional[int], Optional[int], Optional[int]]]] = {}

    def _strip_comments(line: str) -> str:
        # SU2 commonly uses '%' for comments; also tolerate '#'
        line = line.split("%", 1)[0]
        line = line.split("#", 1)[0]
        return line.strip()

    def _parse_keyval(line: str) -> Tuple[str, str]:
        # e.g. "MARKER_TAG= inlet"
        left, right = line.split("=", 1)
        return left.strip(), right.strip()

    current_marker: Optional[str] = None
    remaining_marker_elems: int = 0

    with filepath.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = _strip_comments(raw)
            if not line:
                continue

            # If we are currently consuming marker element lines
            if current_marker is not None and remaining_marker_elems > 0:
                parts = line.split()
                if len(parts) < 3:
                    raise ValueError(
                        f"Malformed marker element line in {filepath} under marker '{current_marker}': {raw!r}"
                    )

                elem = int(parts[0])
                nodes = [int(x) for x in parts[1:]]

                # Keep exactly 3 node columns (p1,p2,p3), pad with NA if needed
                p1 = nodes[0] if len(nodes) > 0 else None
                p2 = nodes[1] if len(nodes) > 1 else None
                p3 = nodes[2] if len(nodes) > 2 else None

                marker_rows.setdefault(current_marker, []).append((elem, p1, p2, p3))

                remaining_marker_elems -= 1
                if remaining_marker_elems == 0:
                    current_marker = None
                continue

            # Otherwise, look for marker headers
            if "=" in line:
                key, val = _parse_keyval(line)

                if key == "MARKER_TAG":
                    current_marker = val
                    marker_rows.setdefault(current_marker, [])
                    # We don't start reading elems until we see MARKER_ELEMS for this tag
                    remaining_marker_elems = 0

                elif key == "MARKER_ELEMS":
                    if current_marker is None:
                        # Some files might list MARKER_ELEMS without immediately preceding MARKER_TAG
                        # This is unusual; treat as format error.
                        raise ValueError(
                            f"Found MARKER_ELEMS before MARKER_TAG in {filepath}: {raw!r}"
                        )
                    remaining_marker_elems = int(val)

                # Ignore other keys (NDIME, NPOIN, NELEM, etc.)
                continue

    # Build DataFrames
    dfs: List[pd.DataFrame] = []
    for marker, rows in marker_rows.items():
        df = pd.DataFrame(rows, columns=["elem", "p1", "p2", "p3"])

        # Use pandas nullable integer dtype so missing p3 is allowed
        df["elem"] = df["elem"].astype("Int64")
        df["p1"] = df["p1"].astype("Int64")
        df["p2"] = df["p2"].astype("Int64")
        df["p3"] = df["p3"].astype("Int64")

        df.name = marker  # store marker tag name on the dataframe
        dfs.append(df)

    return dfs


def su2_markers_dict(filepath: Union[str, Path]) -> Dict[str, pd.DataFrame]:
    """
    Convenience wrapper if you'd rather have a dict keyed by marker tag.
    """
    dfs = read_su2_markers_as_dfs(filepath)
    return {df.name: df for df in dfs}


def get_unique_points(df: pd.DataFrame)-> list:
    cols = [c for c in ("p1", "p2", "p3") if c in df.columns]

    pts = (
        pd.unique(df[cols].to_numpy().ravel())   # unique over p1,p2,p3
    )

    # remove NaNs / <NA> if present, cast to int, sort
    pts = [int(x) for x in pts if pd.notna(x)]
    pts.sort()
    return pts
    

def filter_df_by_point_ids(df: pd.DataFrame, pts: [int], copy_df:bool = False)-> pd.DataFrame: 
    pts_set = set(pts)

    mask = df["PointID"].isin(pts_set)
    out = df.loc[mask]

    return out.copy() if copy_df else out


# =============================================================================
# 
# =============================================================================
filepath = Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\BackPressure_SI_10_tot")
csv_name = "surface_flow_markers.csv"


df = pd.read_csv(filepath / csv_name)


# df_mesh_surfaces = su2_markers_dict(Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\BackPressure_SI_01\Mesh.su2"))
# af_pts = get_unique_points(df_mesh_surfaces["Airfoils"])
# df_af = filter_df_by_point_ids(df_all_data, af_pts, copy_df=True)

# df2 = df2["Airfoils"]



tol = 1e-6

x_LE = 5 * 0.0254
x_TE = 7 * 0.0254
r_hub = 0.2 * 0.0254 
r_duc = 1.825 * 0.0254

# Cut down LE and TE
df = df[(df["x"] > x_LE)] 
df = df[(df["x"] < x_TE)]

# Cut out outter radii
df = df[(df["y"]**2 + df["z"]**2) > r_hub**2 + tol]
df = df[(df["y"]**2 + df["z"]**2) < r_duc**2 - tol]

df.to_csv(filepath / "airfoil_data.csv", index=False)


fig = plt.figure()
ax = fig.add_subplot(projection='3d') # Use "projection='3d'" to create 3D axes

# 3. Plot the data using ax.scatter()
ax.scatter(df['x'], df['y'], df['z'], marker="o")


