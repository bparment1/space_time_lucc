# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
F:\Apps\WinPython-32bit-2.7.6.3\settings\.spyder2\.temp.py
"""
#!/usr/bin/python
#
###################### GENERATE SPATIAL PREDICITONS  ########################
#
#Script to produce spatial predictions with sptial lag/error model.
#
## TODO:
#  - Add paralelization...
#  - get statistics for layers in tif
#
# Authors: Sean McFall, Benoit Parmentier 
# Created on: 11/04/2014
# Updated on: 01/23/2015
# Project: WM Space Beats Time
#
############### LOAD LIBRARY/MODULES USED IN THE SCRIPT #####################

#from pysal import pysal as ps
import pysal as ps

import numpy as np
import os
#import glob

################ NOW FUNCTIONS  ###################

##------------------
# Functions used in the script 
##------------------

def create_dir_and_check_existence(path):
    #Create a new directory
    try:
        os.makedirs(path)
        
    except:
        print "directory already exists"

#######################################################################
######################### BEGIN SCRIPT  ###############################
            
#def main():
    #
########## READ AND PARSE PARAMETERS AND ARGUMENTS ######### 

np.random.seed(100)

ds_name = "EDGY"
y_name = "v1"
x_name = "something else"
out_suffix = "_01232015"

#rootdir = '/home/sean/Documents/Research/space_time_lucc/output/'
#srootdir = '/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/output__predictions_09252014' #on Atlas
#on Benoit's mac
#rootdir = '/Users/benoitparmentier/Google Drive/Space_beats_time/stu/Katrina/output__predictions_09252014' #on 
inDir ="/Users/benoitparmentier/Google Drive/Space_beats_time/output_EDGY_predictions_01232015"
outDir = inDir

#out_dir = "output_data_"+out_suffix
#out_dir = os.path.join(in_dir,out_dir)
#create_dir_and_check_existence(out_dir)
            
#os.chdir(out_dir)        #set working directory


mlerr_dict = {}
# unique id column
unique_id = np.linspace(1, 373, num = 373)

#r_poly_t_155EDGY_predictions_01232015
for date in range(153,156):
    
    #in_file = os.path.join(rootdir, "r_poly_t_" + str(year) + "_predictions_09252014.dbf")
    in_file = os.path.join(inDir, "r_poly_t_" + str(date) + "EDGY_predictions"+ out_suffix +".dbf")
    #r_poly_t_2004_predictions_09252014.dbf
    #katPred = ps.open(rootdir + "r_poly_t_" + str(year) + "_predictions.dbf")
    #r_poly_t_2004_predictions_09252014.dbf
    katPred = ps.open(in_file)
    
    y = np.array([katPred.by_col(y_name)]).T
    
    # array of random numbers 
    x = np.random.randn(len(y),1)
    
    # array of ones
    #x = np.ones_like(y, dtype = np.int)
    #r_poly_t_2004_predictions_09252014
    in_file = os.path.join(inDir, "r_poly_t_" + str(date) + "EDGY_predictions"+ out_suffix +".gal")
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
    out_file = os.path.join(outDir,"mlerrSummary" + str(year) +out_suffix+ ".txt")
    txt_file = open(out_file, "w")
    txt_file.write(mlerr.summary)
    txt_file.close()
    
    # output to array 
    
    # fourth column is the randomly generated numbers
    all_data = np.column_stack((unique_id, mlerr.u, y[:,0], x[:,0], mlerr.predy))
    
    # why does the header insert a '#' symbol onto the front of the string? 
    out_file = os.path.join(outDir,"pysal_res_" + str(year) + out_suffix + ".csv")
    np.savetxt(out_file, all_data, delimiter = ',', 
               header = "unique_ID," + str(year) + "_res,ObsValues,randomVals,predVals")
    
############################### END OF SCRIPT ###########################################

#x_coords = np.genfromtxt(rootdir + "dat_out__predictions.txt", usecols = (1), delimiter = ',', dtype = None)
#y_coords = np.genfromtxt(rootdir + "dat_out__predictions.txt", usecols = (2), delimiter = ',', dtype = None)
#
#x , y = np.loadtxt(rootdir + "dat_out__predictions.txt", delimiter = ',', usecols = (1,2), unpack = True, skiprows = 1)
#

# output dictionary to array file
