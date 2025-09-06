# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 08:20:29 2023

@author: cycon
"""

'''
File Description
I want this module to generate a set of vertices (and also output faces) 
for a rotor based on certain inputs:
    - Root airfoil
        - With chord width cr
    - Tip Airfoil
        - With chord width ct
    - Twist Angle
    - Root AoA
    - Hub Diameter
    - Hub height
    - Blade height
    - Number of blades
Plan:
    Gather inputs and import airfoils as Point Clouds.
    Must ensure that the points are in the right order (start at LE, end at LE)
        - These will be 2d in XY plane. I plan on having the z-axis align with hub center
    
    
    In order to avoid twist-tightening (when twisting between two profiles results
                                        in a shape with too thin of a cross-section
                                        due to the twist), I think that I can use 
        an iterative method of producing a large amount of sub-aerofoils in between
        the root and tip. This would work by interpolating between the two airfoils 
        based on its relative radius to blade-height position. I can do this between
        two airfoils with non-consistent x-spacing by using an average value at a particular 
        dx. I may want to apply this to each airfoil prior to beginning to ensure that 
        each airfoil is described by the same number of points that correlate to the same
        x-values (ie point 50 on two different airfoils describes the top-surface at x/l = 0.5)
        This will also allow me to auto-generate faces by following the curve all the way around.
        in between each sub-airfoil.
        
'''