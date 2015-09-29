# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
F:\Apps\WinPython-32bit-2.7.6.3\settings\.spyder2\.temp.py
"""

#from pysal import pysal as ps
import pysal as ps

import numpy as np
import os
#import glob

np.random.seed(100)

ds_name = "katrina2004"
y_name = "v1"
x_name = "something else"
out_suffix = "_09252014"

#rootdir = '/home/sean/Documents/Research/space_time_lucc/output/'
rootdir = '/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/output__predictions_09252014'

mlerr_dict = {}
# unique id column
unique_id = np.linspace(1, 373, num = 373)
for year in range(2001,2013):
    
    #in_file = os.path.join(rootdir, "r_poly_t_" + str(year) + "_predictions_09252014.dbf")
    in_file = os.path.join(rootdir, "r_poly_t_" + str(year) + "_predictions"+ out_suffix +".dbf")
    #katPred = ps.open(rootdir + "r_poly_t_" + str(year) + "_predictions.dbf")
    #r_poly_t_2004_predictions_09252014.dbf
    katPred = ps.open(in_file)
    
    y = np.array([katPred.by_col(y_name)]).T
    
    # array of random numbers 
    x = np.random.randn(len(y),1)
    
    # array of ones
    #x = np.ones_like(y, dtype = np.int)
    #r_poly_t_2004_predictions_09252014
    in_file = os.path.join(rootdir, "r_poly_t_" + str(year) + "_predictions"+ out_suffix +".gal")
    ww = ps.open(in_file, 'r')

    #ww = ps.open(rootdir + "r_poly_t_" + str(year) + "_predictions" + out_suffix +".gal", 'r')
    w = ww.read()
    ww.close()
    w_name = str(year) + "pred.gal"
    w.transform = 'r'
    
    print year
    mlerr = ps.spreg.ml_error.ML_Error(y, x, w,method="ord")
    #mlerr = ps.spreg.ml_error.ML_Error(y, x, w)
    
    mlerr_dict[year] = mlerr
    
    # output summary to individual text file
    txt_file = open("mlerrSummary" + str(year) + ".txt", "w")
    txt_file.write(mlerr.summary)
    txt_file.close()
    
    # output to array 
    
    # fourth column is the randomly generated numbers
    all_data = np.column_stack((unique_id, mlerr.u, y[:,0], x[:,0], mlerr.predy))
    
    # why does the header insert a '#' symbol onto the front of the string? 
    np.savetxt("pysal_res_" + str(year) + ".csv", all_data, delimiter = ',', 
               header = "unique_ID," + str(year) + "_res,ObsValues,randomVals,predVals")
    

#x_coords = np.genfromtxt(rootdir + "dat_out__predictions.txt", usecols = (1), delimiter = ',', dtype = None)
#y_coords = np.genfromtxt(rootdir + "dat_out__predictions.txt", usecols = (2), delimiter = ',', dtype = None)
#
#x , y = np.loadtxt(rootdir + "dat_out__predictions.txt", delimiter = ',', usecols = (1,2), unpack = True, skiprows = 1)
#

# output dictionary to array file
