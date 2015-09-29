# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
F:\Apps\WinPython-32bit-2.7.6.3\settings\.spyder2\.temp.py
"""

import pysal as ps
import numpy as np


#
# below is the example provided for the spreg.ml_error model
# NOTE: it should work wherever given the relative file paths
#

#db = ps.open(ps.examples.get_path("south.dbf"), 'r')
#ds_name = "south.dbf"
#y_name = "HR90"
#y = np.array(db.by_col(y_name))
#y.shape = (len(y), 1)
#x_names = ["RD90","PS90","UE90","DV90"]
#x = np.array([db.by_col(var) for var in x_names]).T
#ww = ps.open(ps.examples.get_path("south_q.gal"))
#w = ww.read()
#ww.close()
#w_name = "south_q.gal"
#w.transform = 'r'
#mlerr = ps.spreg.ml_error.ML_Error(y,x,w,name_y=y_name,name_x=x_names,name_w=w_name,name_ds=ds_name)

#
# below is my take on the code with the SBT data
#

# converted the txt file to csv, pysal seems to like the file type more
katPred = ps.open("F:\Research\spaceTime\output\slm_2014-10-30\dat_out__predictions2014-10-30.csv", 'r')

ds_name = "katrina2004"
y_name = "LSPOP2004"
x_name = y_name + "x"

# y = np.array(katPred.by_col(y_name))

# given that we do not want a covariate, I made the DV and IV the same
# must be a nx1 array to work
x = np.array([katPred.by_col(y_name)]).T
y = x

# does not like number in the file name, though it's ok with the directory name???
# I am guessing this is a windows problem
ww = ps.open("F:\Research\pred.gal", 'r')
w = ww.read()
ww.close()
w_name = "2004pred.gal"
# not sure what the below line of code does, but it was in the example...
w.transform = 'r'

mlerr = ps.spreg.ml_error.ML_Error(y, x, w, name_y = y_name, name_x = x_name, name_w = w_name, name_ds = ds_name)
