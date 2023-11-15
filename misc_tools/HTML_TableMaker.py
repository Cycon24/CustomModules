# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:40:43 2023

@author: cycon
"""
import numpy as np
from IPython.display import Markdown as md

'''
This module will have a function that returns the string containing HTML
that will be able to print markdown talbes with values from python
utilizing .format()
'''


# Needs some revision. Needs to be able to handle (n,) arrays
def checkSize(toCheck, nRow, nCol):
    matrix = np.array(toCheck)
    # print(np.shape(matrix))
    if np.shape(matrix) == (nRow,nCol):
        return True
    else:
        return False

def formatCorrector(inFormatter, size):
    defaultForm = '{}'
    [rows, cols] = size
    outFormatter = np.empty((rows,cols),dtype="object")
    
    if inFormatter == None:
        for i in range(0,rows):
            for j in range(0,cols):
                # Since no format was given, print the values as default
                outFormatter[i,j] = defaultForm
        return 'None', outFormatter #________________
            
    # Check if a single value was inputted -> fill whole array w value
    elif type(inFormatter) != np.ndarray and type(inFormatter) != list:
        for i in range(0,rows):
            for j in range(0,cols):
                # Since a single format was given, use the value for all of them
                outFormatter[i,j] = inFormatter
        return 'Single Value', outFormatter #________________
    
    # The format is an array
    else:
        if len(np.shape(inFormatter)) == 1: # So inf is in shape (nrows,) to change to (1,ncols=nrows) since we view [x1, x2, ... xn] as a 1row vec
            inFormatter = np.reshape(inFormatter, (1,np.shape(inFormatter)[0]))
        else:
            fr,fc = np.shape(inFormatter)
            inFormatter = np.reshape(inFormatter, (fr,fc))
            
        
        # Check the size
        if checkSize(inFormatter, rows, cols):
            # Same size
            outFormatter = inFormatter
            return 'Array - Same size', outFormatter
            
        else:
            # Different size, need to check what axis isnt same size
            frows, fcols = np.shape(inFormatter)
            if frows == 1 and fcols < cols:
                # Going to repeat row values down column for rows we have
                for j in range(0,fcols):
                    outFormatter[0 ,j] = inFormatter[0,j]
                
                # First fill rest of cols with empty
                for j in range(fcols, cols):
                        outFormatter[0,j] = defaultForm
                
                # Fill rest of rows with copy of 1st
                for i in range(1, rows):
                    outFormatter[i,:] = outFormatter[0,:]
                
                return 'Array 1xSmaller', outFormatter
                
            elif fcols == 1 and frows < rows:
                # Going to repeat col values across rows for cols we have
                for i in range(0,frows):
                    outFormatter[i ,0] = inFormatter[i,0]
                
                
                # First fill rest of rows with empty
                for i in range(frows, rows):
                        outFormatter[i,0] = defaultForm
                
                # Fill rest of rows with copy of 1st
                for j in range(1, cols):
                    outFormatter[:,j] = outFormatter[:,0]
                
                return 'Array Smallerx1', outFormatter

            elif frows < rows and fcols < cols:
                # fill outform with inform for what we do have
                # for j in range(0,fcols):
                #     for i in range(1, rows):
                #         outFormatter[i,j] = outFormatter[i,j]
                    
                    
                outFormatter[0:frows,0:fcols] = inFormatter
                # Fill rest of outForm with blank formats
                for i in range(0,rows):
                    for j in range(0, cols):
                        if i >= frows or j >= fcols:
                            outFormatter[i,j] = defaultForm
                return 'Array - smallerxsmaller', outFormatter
           
            elif frows == 1 and fcols == cols:
                # Going to repeat row values down column for rows we have
                for j in range(0,cols):
                    outFormatter[0 ,j] = inFormatter[0,j]
                
                # Fill rest of rows with copy of 1st
                for i in range(1, rows):
                    outFormatter[i,:] = outFormatter[0,:]
                
                return 'Array 1xSame', outFormatter
                
            elif fcols == 1 and frows == rows:
                # Going to repeat col values across rows for cols we have
                for i in range(0,rows):
                    outFormatter[i ,0] = inFormatter[i,0]
                
                # Fill rest of rows with copy of 1st
                for j in range(1, cols):
                    outFormatter[:,j] = outFormatter[:,0]
                
                return 'Array Samex1', outFormatter
            
            else:
                print('ERROR: Formatting cases not satisfied')
                return None, None

'''
Make Table will return a markdown object (Or String) that when ran in 
a jupyter notebooke will generate a pretty HTML table. colHeaders and colRows
are single dimension array of strings (or numbers) that will be the headers. 
data can be a 2D array in which the rows/columns corresponds to the col and row
headers.
Currently does not support row headers being none
'''
def makeTable(colHeaders, rowHeaders, data, formatter = None, asString=False):
    tableStr = '<table>'
    data = np.array(data)
    
    if rowHeaders == None:
        noRowHeader = True
        rows=len(data[:,0])
    else:
        noRowHeader = False
        rows = len(rowHeaders)
    
    cols = len(colHeaders)
    
    
    if len(np.shape(data)) == 1: # So data is in shape (nrows,) to change to (1,ncols=nrows) since we view [x1, x2, ... xn] as a 1row vec
        data = np.reshape(data, (1,np.shape(data)[0]))
        
    # Check data size
    if not checkSize(data, rows, cols):
        print('Error: Data for table cannot fit shape of row and column headers.')
        print('\tCannot fit {} data into {}x{} table'.format(np.shape(data),rows,cols))
        return None
    
    # Check formatter
    frmStrCheck, formatter = formatCorrector(formatter, [rows, cols])
    
    # Add Table Title
    
    # Add Col Headers
    tableStr += '<tr><th><th>' if not noRowHeader else '<tr>'
    for i in range(cols):
        tableStr += '<th> {} <th>'.format(colHeaders[i])
    tableStr += '<tr>'
    
    # Add Rows
    for i in range(rows):
        tableStr += '<tr>'
        
        # Row Header
        if not noRowHeader:
            tableStr += '<th> {} <th>'.format(rowHeaders[i])
        
        # Add table data? May need to pass data in
        for j in range(cols):
            formattedData = formatter[i,j].format(data[i,j])
            tableStr+= '<td>{}<td>'.format(formattedData)
            
        tableStr += '<tr>'
        
    tableStr += '<table>'
    return md(tableStr) if not asString else tableStr
        
'''    
# Example Table
### Results Table\
<table>\
   <tr> \
       <th colspan="7">Values Averaged over {10} Rounds <th>\
       <th><th>\
       <th><th>\
   <tr>\
   <tr>\
       <th><th>\
       <th> Gauss Seidel <th>\
       <th> Conj. Gradient <th>\
       <th> Gauss Pivot <th>\
       <th> LU Pivot <th>\
       <th> Cramer <th>\
   <tr>\
   <tr>\
       <th> Ave Compute Times <th>\
       <td> {0:8.6f} ms <td>\
       <td> {1:8.6f} ms <td>\
       <td> {2:8.6f} ms <td>\
       <td> {3:8.6f} ms <td> \
       <td> {4:8.6f} ms <td> \
    <tr>\
    <tr>\
       <th> Averaged Residual <th>\
       <td> {5:10.6E} <td>\
       <td> {6:10.6E} <td>\
       <td> {7:10.6E} <td>\
       <td> {8:10.6E} <td>\
       <td> {9:10.6E} <td>\
    <tr>\
<table>'
'''


    

if __name__ == "__main__":
    # formstr = "{:.3f}"
    # teststr = "Hello this is a test {}"
    # print(teststr)
    # teststr = teststr.format(formstr)
    # print(teststr.format(1.54))
    # print(checkSize([[3,1,4]], 1, 3))
    # print(type(np.zeros(1)))
    # print(np.reshape(np.zeros((5,1)), (1,5)))
    
    # first = np.zeros((4,6))
    # second = np.ones((7,9))
    # second[0:4,0:6] = first
    # print(second)
    # for i in range(0,7):
    #     for j in range(0,9):
    #         if i >= 4 or j >= 6:
    #             second[i,j] = 2
    # print(second)
    # a = np.zeros(4)
    # print(np.shape(a),a)
    # a = np.reshape(a, (1,np.shape(a)[0]))
    # print(np.shape(a),a)
    # print(type(['{:2.f}','{:.3f}'])==list)
    
    # ---- Format Case Testing -----
    # None type
    # strreturn, arrreturn = formatCorrector(None, [3,3])
    # print( strreturn,'\n', arrreturn)
    
    # # 1 row by smaller cols (2,)->(1,2)
    # strreturn, arrreturn = formatCorrector(['{:2.f}','{:.3f}'], [4,3])
    # print( strreturn,'\n', arrreturn,'\n')
    
    # # 1 row by smaller cols (1,2)
    # strreturn, arrreturn = formatCorrector([['{:.2f}', '{:.3f}']], [20,3])
    # print( strreturn,'\n', arrreturn,'\n')
    
    # # 1 col by smaller rows (2,1)
    # strreturn, arrreturn = formatCorrector([['{:2.f}'],['{:.3f}']], [4,3])
    # print( strreturn,'\n', arrreturn,'\n')
    
    # # smaller by smaller (2,2)
    # strreturn, arrreturn = formatCorrector([['{:2.f}','{:.3f}'],['{:4.f}','{:.5f}']], [4,3])
    # print( strreturn,'\n', arrreturn,'\n')
    # # print(None==np.shape(a)[1])
    
    # #Same size (3,3)
    # strreturn, arrreturn = formatCorrector([['{:1.f}','{:2.f}','{:.3f}'],['{:4.f}','{:.5f}','{:.6f}'],['{:7.f}','{:.8f}','{:.9f}']], [3,3])
    # print( strreturn,'\n', arrreturn,'\n')
    # # print(None==np.shape(a)[1])
    
    
    # #1 x Same size (1,3)
    # strreturn, arrreturn = formatCorrector(['{:5.0f}','{:5.1f}','{:5.4f}'], [3,3])
    # print( strreturn,'\n', arrreturn,'\n')
    
    
    # ---------------------
    # ----- Table Testing -----
    # colHeaders = ['First C', 'Second C', 'Third C']
    # rowHeaders = ['First R', 'Second R', 'Third R']

    # data = np.zeros((3,3))
    # for i in range(0,3):
    #     for j in range(0,3):
    #         data[i,j] = i + j
    # tabStr = makeTable(colHeaders, rowHeaders, data)
    # print(tabStr)
    
    # ------ No Row Header Testing -----
    # myStr = '1'
    # myStr += '2' if 1==2 else '3'
    # print(myStr)
    
    # headers = ['alpha','calc1','calc2']
    
    # alpha = np.arange(0,20,1)
    # calc1 = np.ones(np.size(alpha))
    # calc2 = np.ones(np.size(alpha))
    # calc1 = np.multiply(calc1,alpha)*5
    # calc2 = np.multiply(calc2,alpha)*np.pi
    

    # data = np.zeros((np.size(alpha),3))
    # data[:,0] = alpha
    # data[:,1] = calc1
    # data[:,2] = calc2
    
    
    # tabStr = makeTable(headers, None, data, None, True)
    # print(tabStr)
    
    print('END')