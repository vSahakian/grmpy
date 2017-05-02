## Just make Vs30 plots for now...


import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr


#fig dir:
png_dir = '/Users/vsahakian/anza/models/residuals/v3_vs30_compare/figs/'
pdf_dir = '/Users/vsahakian/anza/models/residuals/v3_vs30_compare/figs/pdf/' 

# Model 1, 6 coefficients, iterative ME with a4 and a5 fixed
me6cff_iter_path = '/Users/vsahakian/anza/models/residuals/mixedregr_v3anza2013_pga_6coeff_a4_-1.73_a5_-0.01_Mc_8.5_res4/mixedcoeff_v3anza2013_pga_vs30_6coeff_pga__ncoeff6_Mc_8.5_VR_99.3_a4_-1.73_a5_-0.01_robj.pckl'

modfile = open(me6cff_iter_path,'r')
me6cff_iter = pickle.load(modfile)
modfile.close()

#######################
# Model 2, 5 coefficients, iterative ME with a4 and a5 fixed
me5cff_iter_path = '/Users/vsahakian/anza/models/residuals/mixedregr_v3anza2013_pga_5coeff_a4_-1.72_a5_-0.01_Mc_8.5_res4_noVs30/mixedcoeff_v3anza2013_pga_noVs30_5coeff_pga__ncoeff5_Mc_8.5_VR_99.4_a4_-1.72_a5_-0.01_robj.pckl'

modfile = open(me5cff_iter_path,'r')
me5cff_iter = pickle.load(modfile)
modfile.close()


###########################################
# Unique stations:
usta,usta_ind = np.unique(me6cff_iter.stnum,return_index=True)



###########################################
# Sites vs. vs30....

vs30fig = plt.figure(figsize=(10,8))
vs30png = png_dir + 'sites_vs30.png'
vs30pdf = pdf_dir + 'sites_vs30.pdf'

# Get correlation values
me6corr,me6r=pearsonr(me6cff_iter.vs30[usta_ind],np.array(me6cff_iter.site_terms)[usta_ind])
me5corr,me5r=pearsonr(me5cff_iter.vs30[usta_ind],np.array(me5cff_iter.site_terms)[usta_ind])

# sort by vs30:
vs30sort = np.argsort(me6cff_iter.vs30[usta_ind])

# Get difference to plot:
me_site_diff_unsort = np.array(me6cff_iter.site_terms)[usta_ind] - np.array(me5cff_iter.site_terms)[usta_ind]
me_site_diff = me_site_diff_unsort[vs30sort]

# Plot..
plt.scatter(me6cff_iter.vs30[usta_ind],np.array(me6cff_iter.site_terms)[usta_ind],s=25,edgecolor='#5d9fba',facecolor='none',marker='D',linewidth=2,label='Mixed iterative, with Vs30, coeff = %.2f, r = %.2f' % (me6corr,me6r)) 
plt.scatter(me5cff_iter.vs30[usta_ind],np.array(me5cff_iter.site_terms)[usta_ind],s=25,edgecolor='#b57955',facecolor='none',marker='o',linewidth=2,label='Mixed iterative, without Vs30, coeff = %.2f, r = %.2f' % (me5corr, me5r)) 

plt.plot(me6cff_iter.vs30[usta_ind][vs30sort], me_site_diff,linestyle='--',linewidth=2,color='#4f4f4f',label='Residual (vs30 - no vs30)')

plt.ylim([-2,3.5])

plt.xlabel('Vs30 (m/s)')
plt.ylabel('Site Term (ln residual)')
plt.title('Site term vs. Vs30')

plt.legend(loc=2)

plt.savefig(vs30png)
plt.savefig(vs30pdf)