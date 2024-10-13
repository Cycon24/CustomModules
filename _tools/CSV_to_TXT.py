# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:55:51 2023

@author: cycon
"""

import csv

def convert(inFileName):
    csv_file = inFileName + '.csv'
    txt_file = inFileName + '.txt'
    with open(txt_file, "w") as my_output_file:
        with open(csv_file, "r") as my_input_file:
            [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()
        


# C:\Users\cycon\OneDrive - University of Cincinnati\Spring 2023\Aerodynamics\Airfoils
if __name__=='__main__':
    fileLoc =  'C:\\Users\\cycon\\OneDrive - University of Cincinnati\\Spring \
2023\\Aerodynamics\\Assignments\\Team Project\\AirfoilCurves\\'
    Profiles = ['NACA 2706','NACA 2709']
    
    for p in Profiles:
        convert(fileLoc+p)
    print('Complete')
        