# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 14:23:40 2025

@author: BriceM

This code will be interacting with the previously refined mesh generation code,
SU2 runner, and the configuration files. In order to do this, we need interfaces
for each that includes the parameters we will use as inputs. We want to develop a 
consistent naming system for results and file management, as we will be generating
a large number of files for each sweep. Additionally, we want to maintain a log
file that tracks the overall sweep progression. Below we will define the parameters
necessary for each code.

General Flow (used by mesh and config):
- Mi
- Pi
- Ti
- Gas:
    - gamma
    - R
    - mu 
- Other calculated flow vars (total properties, vel, Re, etc)

Mesh
- Pipe L Upstream
- Pipe L Downstream
- Crossection Radius
- Number of Blades
- Chord Length
- Blade shape
- Blade AoA
- Boundary Layer Parameters
- General Flow vars

Config File
- Iterations
- General Flow vars
- File names
- File output names

SU2 Runner
- File names
- File locations/folders



"""

