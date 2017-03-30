## Plot the differences between models...
# VJS 3/2017

import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt


## Read in models: 

# Model 1: 6 coefficients, single-mean
sm6cff_path = '/Users/vsahakian/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta_res4_pga_vs30_6coeff/v2anza2013_Mc8.5_pgrid_5sta_res4_pga_vs30_6coeff_robj.pckl'

modfile = open(sm6cff_path,'r')
sm6cff = pickle.load(modfile)
modfile.close()

# Model 2: 6 coefficients, ME all
me6cff_path = '/Users/vsahakian/anza/models/residuals/mixedregr_v2anza2013_pga_6coeff_Mc_8.5_res4/mixedcoeff_v2anza2013_pga_vs30_6coeff_pga__ncoeff6_Mc_8.5_VR_99.9_robj.pckl'

modfile = open(me6cff_path,'r')
me6cff = pickle.load(modfile)
modfile.close()

# Model 3, 6 coefficients, iterative ME intercept
me6cff_iter_path = '/Users/vsahakian/anza/models/residuals/mixedregr_v2anza2013_pga_6coeff_a2_0.73_a3_-0.14_a4_-1.68_a5_-0.01_a6_1.03_Mc_8.5_res4/mixedcoeff_v2anza2013_pga_vs30_6coeff_pga__ncoeff6_Mc_8.5_VR_81.4_a2_0.73_a3_-0.14_a4_-1.68_a5_-0.01_a6_1.03_robj.pckl'

modfile = open(me6cff_iter_path,'r')
me6cff_iter = pickle.load(modfile)
modfile.close()

##################

# Model 4: 5 coefficients, single-mean
sm5cff_path = '/Users/vsahakian/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta_res4_pga_noVs30_5coeff/v2anza2013_Mc8.5_pgrid_5sta_res4_pga_noVs30_5coeff_robj.pckl'

modfile = open(sm5cff_path,'r')
sm5cff = pickle.load(modfile)
modfile.close()

# Model 5: 5 coefficients, ME all
me5cff_path = '/Users/vsahakian/anza/models/residuals/mixedregr_v2anza2013_pga_5coeff_Mc_8.5_res4_noVs30/mixedcoeff_v2anza2013_pga_noVs30_5coeff_pga__ncoeff5_Mc_8.5_VR_99.9_robj.pckl'

modfile = open(me5cff_path,'r')
me5cff = pickle.load(modfile)
modfile.close()

# Model 6, 5 coefficients, iterative ME intercept
me5cff_iter_path = '/Users/vsahakian/anza/models/residuals/mixedregr_v2anza2013_pga_5coeff_a2_0.31_a3_-0.17_a4_-1.72_a5_-0.01_Mc_8.5_res4_noVs30/mixedcoeff_v2anza2013_pga_noVs30_5coeff_pga__ncoeff5_Mc_8.5_VR_97.8_a2_0.31_a3_-0.17_a4_-1.72_a5_-0.01_robj.pckl'

modfile = open(me5cff_iter_path,'r')
me5cff_iter = pickle.load(modfile)
modfile.close()


#######################
## Coefficient plots




## Std plots

sm6cff_std = np.array([np.std(sm6cff.total_residual),sm6cff.E_std[0],sm6cff.site_std,sm6cff.path_std])



## Diff. between event term plots




## Diff. between site term plots




## Diff. between path term plots

