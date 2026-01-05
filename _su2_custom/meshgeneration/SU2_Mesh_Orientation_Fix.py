# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 17:30:44 2025

@author: BriceM (Chat CPT provided)
"""

import numpy as np

def read_su2_mesh(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    return lines

def write_su2_mesh(filename, lines):
    with open(filename, "w") as f:
        f.writelines(lines)

def parse_points(lines):
    points = []
    npoint = 0
    for i, line in enumerate(lines):
        if line.startswith("NPOIN="):
            npoint = int(line.split("=")[1].strip())
            for j in range(i + 1, i + 1 + npoint):
                parts = lines[j].split()
                coords = list(map(float, parts[:2]))  # only 2D coords
                points.append(coords)
            break
    return np.array(points)

def fix_tri_orientation(lines, area_tol=1e-6):
    points = parse_points(lines)
    fixed_lines = []
    WITHIN_ELEMENTS = False
    num_fixed_elems = 0
    num_elems_checked = 0
    num_elems_removed = 0
    for line in lines:
        if line.startswith("NPOIN="):
            WITHIN_ELEMENTS = False
        
        if line.startswith("NELEM="):
            WITHIN_ELEMENTS = True
        
        if line.strip().startswith("5") and WITHIN_ELEMENTS:
            num_elems_checked += 1
            parts = line.split()
            n1, n2, n3 = map(int, parts[1:4])
            p1, p2, p3 = points[n1], points[n2], points[n3]

            # Compute signed area
            area = 0.5 * ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]))
            
            # Negative area, flip points 
            if area < 0:
                # Flip orientation (swap last 2 nodes)
                n2, n3 = n3, n2
                num_fixed_elems += 1
                
            # Area too small
            # if abs(area) < area_tol:
            #     # if its too small, move to next iteration without writing
            #     num_elems_removed+=1
            #     continue  

            fixed_line = f"5 {n1} {n2} {n3} {parts[4]}\n"
            fixed_lines.append(fixed_line)
        else:
            fixed_lines.append(line)
    print(f'Checked {num_elems_checked} triangles...')
    print(f'Fixed {num_fixed_elems} triangles...')
    print(f'Removed {num_elems_removed} triangles...')
    
    return fixed_lines

if __name__ == "__main__":
    infile = "mesh_in.su2"
    outfile = "mesh_fixed.su2"

    lines = read_su2_mesh(infile)
    fixed_lines = fix_tri_orientation(lines)
    write_su2_mesh(outfile, fixed_lines)
    print(f"Fixed mesh written to {outfile}")
