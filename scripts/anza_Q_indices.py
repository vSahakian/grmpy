#############  Script to compute Q metrics  ###########
## VJS 2/2018

import dread
import cdefs as cdf
import cPickle as pickle
from numpy import meshgrid,zeros,where,array
import run_res_analysis as runra
import res_analysis as ra
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

#######################################
##############  Paths #################
#######################################

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'


## Path to Qp and Qs
hauksson_model_path_Qp = HOME + '/anza/data/vm/Hauksson/Qp.pckl'
hauksson_model_path_Qs = HOME + '/anza/data/vm/Hauksson/Qs.pckl'


