#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:25:53 2020

@author: vjs
"""

## For Amanda - get GMM estimations for So Cal small M

import cdefs
import dill
import pickle as pickle
import os

p2_modelpath = '/Users/vjs/anza/models/pckl/v3anza2013_pga_noVs30_5coeff_a4_-1.2/mixedregr_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2.pckl'

output_dir = '/Users/vjs/anza/models/pckl/v3anza2013_pga_noVs30_5coeff_a4_-1.2/'

def convertp2top3(old_pkl,output_directory):
    """
    Convert a Python 2 pickle to Python 3
    """
    # Make a name for the new pickle
    new_pkl = output_directory + os.path.splitext(os.path.basename(old_pkl))[0]+"_p3.pkl"

    # Convert Python 2 "ObjectType" to Python 3 object
    dill._dill._reverse_typemap["ObjectType"] = object

    # Open the pickle using latin1 encoding
    with open(old_pkl, "rb") as f:
        loaded = pickle.load(f, encoding="latin1")

    # Re-save as Python 3 pickle
    with open(new_pkl, "wb") as outfile:
        pickle.dump(loaded, outfile)

    print('Wrote to file ' + new_pkl)

### Convert it to a python 3 pickle:
    
convertp2top3(p2_modelpath,output_dir)