#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 16:10:00 2020

@author: vjs
"""


import numpy as np
import gmpe as gm
import pickle as pickle
import matplotlib.pyplot as plt

## Path to the pickle file with model parameters
pickle_path = '/Users/vjs/anza/models/pckl/v3anza2013_pga_noVs30_5coeff_a4_-1.2/mixedregr_v3anza2013_pga_noVs30_5coeff_a4_-1.2_pga__ncoeff5_Mc_8.5_VR_99.6_a4_-1.2_p3.pkl'

## Mag to compute at:
mw = 1.0

## Number of model coefficients
numcoeff = 5

#Use magnitude-dependent fictitious depth/finite fault dimension factor?####
#no == 0, yes == 1
mdep_ffdf = 0

##Centering Magnitude for x**2 term
Mc = 8.5

sdist = np.logspace(np.log10(1),np.log10(250),50)

with open(pickle_path,"rb") as f:
    model = pickle.load(f)

#mw_model,estimated_pga = gm.compute_model_fixeddist(model.m,[1,1.1],sdist,Mc,mdep_ffdf,ncoeff=numcoeff)
estimated_pga = gm.compute_model_fixedmag(model.m,mw,sdist,Mc,mdep_ffdf,ncoeff=5)
#gm.compute_model(model.m,rng,mw,r,ffdf,vs30,Mc,vref,mdep_ffdf,ncoeff)

plt.figure()
plt.plot(sdist,np.exp(estimated_pga))
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distance (km)')
plt.ylabel('PGA (g)')