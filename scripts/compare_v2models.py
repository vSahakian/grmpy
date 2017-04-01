## Plot the differences between models...
# VJS 3/2017

import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt


#fig dir:
png_dir = '/Users/vsahakian/anza/models/residuals/sm_me_iterative_compare/figs/'
pdf_dir = '/Users/vsahakian/anza/models/residuals/sm_me_iterative_compare/figs/pdf/' 

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

sm6cff_std = np.array([np.std(sm6cff.total_residual),sm6cff.E_std,sm6cff.site_std,sm6cff.path_std])



## Diff. between event term plots




#### Diff. between site term plots

## First, plot site terms + errorbars on one plot

# Get unique sites to plot:
usta,usta_ind = np.unique(sm6cff.stnum,return_index=True)

####
# Inititate first:
allsitefig = plt.figure(figsize=(25,8))
allsitefig_png = png_dir + 'all_sites.png'
allsitefig_pdf = pdf_dir + 'all_sites.pdf'


# sm6cff
plt.errorbar(sm6cff.stnum[usta_ind],np.array(sm6cff.site_terms)[usta_ind],yerr=np.array(sm6cff.site_stderr)[usta_ind],fmt='o',color='#687da0',label='Single-mean 6 coefficients')

# me6cff
plt.errorbar(me6cff.stnum[usta_ind],np.array(me6cff.site_terms)[usta_ind],yerr=np.array(me6cff.site_stderr)[usta_ind],fmt='^',markersize=8,color='#4b7763',elinewidth=2,capsize=4,capthick=1.5,label='Mixed 6 coefficients')

# me6cff iterative
plt.errorbar(me6cff_iter.stnum[usta_ind],np.array(me6cff_iter.site_terms)[usta_ind],yerr=np.array(me6cff_iter.site_stderr)[usta_ind],fmt='*',markersize=8,color='#dbb36d',elinewidth=2,capsize=4,capthick=1.5,label='Iterative-Mixed 6 coefficients')

# sm5cff
plt.errorbar(sm5cff.stnum[usta_ind],np.array(sm5cff.site_terms)[usta_ind],yerr=np.array(sm5cff.site_stderr)[usta_ind],fmt='h',color='#bc0303',label='Single-mean 5 coefficients')

# me5cff
plt.errorbar(me5cff.stnum[usta_ind],np.array(me5cff.site_terms)[usta_ind],yerr=np.array(me5cff.site_stderr)[usta_ind],fmt='D',markersize=8,color='#796899',elinewidth=2,capsize=4,capthick=1.5,label='Mixed 5 coefficients')

#m e5cff iterative
plt.errorbar(me5cff_iter.stnum[usta_ind],np.array(me5cff_iter.site_terms)[usta_ind],yerr=np.array(me5cff_iter.site_stderr)[usta_ind],fmt='p',markersize=8,color='#d17c21',elinewidth=2,capsize=4,capthick=1.5,label='Iterative-Mixed 5 coefficients')


# Station name text:
for stationi in range(len(usta)):
    statext = sm6cff.sta[usta_ind][stationi]
    stx = sm6cff.stnum[usta_ind][stationi]
    sty = -2.5
    rotangle = 90
    
    plt.text(stx,sty,statext,rotation=rotangle,ha='center',va='bottom')

plt.legend()
plt.title('Site vs. Site Terms, all runs')

plt.savefig(allsitefig_png)
plt.savefig(allsitefig_pdf)


###########################################################3
# Inititate first:
c6sitefig = plt.figure(figsize=(25,8))
c6sitefig_png = png_dir + 'cff6_sites.png'
c6sitefig_pdf = pdf_dir + 'cff6_sites.pdf'


# sm6cff
plt.errorbar(sm6cff.stnum[usta_ind],np.array(sm6cff.site_terms)[usta_ind],yerr=np.array(sm6cff.site_stderr)[usta_ind],fmt='o',color='#687da0',label='Single-mean 6 coefficients')

# me6cff
plt.errorbar(me6cff.stnum[usta_ind],np.array(me6cff.site_terms)[usta_ind],yerr=np.array(me6cff.site_stderr)[usta_ind],fmt='^',markersize=8,color='#4b7763',elinewidth=2,capsize=4,capthick=1.5,label='Mixed 6 coefficients')

# me6cff iterative
plt.errorbar(me6cff_iter.stnum[usta_ind],np.array(me6cff_iter.site_terms)[usta_ind],yerr=np.array(me6cff_iter.site_stderr)[usta_ind],fmt='*',markersize=8,color='#dbb36d',elinewidth=2,capsize=4,capthick=1.5,label='Iterative-Mixed 6 coefficients')


# Station name text:
for stationi in range(len(usta)):
    statext = sm6cff.sta[usta_ind][stationi]
    stx = sm6cff.stnum[usta_ind][stationi]
    sty = -2.5
    rotangle = 90
    
    plt.text(stx,sty,statext,rotation=rotangle,ha='center',va='bottom')


plt.legend()
plt.title('Site vs. Site terms, 6 coefficients/with Vs30 term only')

plt.savefig(c6sitefig_png)
plt.savefig(c6sitefig_pdf)



###########################################################3
# Inititate first:
c5sitefig = plt.figure(figsize=(25,8))
c5sitefig_png = png_dir + 'cff5_sites.png'
c5sitefig_pdf = pdf_dir + 'cff5_sites.pdf'


# sm5cff
plt.errorbar(sm5cff.stnum[usta_ind],np.array(sm5cff.site_terms)[usta_ind],yerr=np.array(sm5cff.site_stderr)[usta_ind],fmt='h',color='#bc0303',label='Single-mean 5 coefficients')

# me5cff
plt.errorbar(me5cff.stnum[usta_ind],np.array(me5cff.site_terms)[usta_ind],yerr=np.array(me5cff.site_stderr)[usta_ind],fmt='D',markersize=8,color='#796899',elinewidth=2,capsize=4,capthick=1.5,label='Mixed 5 coefficients')

#m e5cff iterative
plt.errorbar(me5cff_iter.stnum[usta_ind],np.array(me5cff_iter.site_terms)[usta_ind],yerr=np.array(me5cff_iter.site_stderr)[usta_ind],fmt='p',markersize=8,color='#d17c21',elinewidth=2,capsize=4,capthick=1.5,label='Iterative-Mixed 5 coefficients')


# Station name text:
for stationi in range(len(usta)):
    statext = sm5cff.sta[usta_ind][stationi]
    stx = sm5cff.stnum[usta_ind][stationi]
    sty = -2.5
    rotangle = 90
    
    plt.text(stx,sty,statext,rotation=rotangle,ha='center',va='bottom')


plt.legend()
plt.title('Site vs. Site terms, 5 coefficients/no Vs30 only')

plt.savefig(c5sitefig_png)
plt.savefig(c5sitefig_pdf)


######################
# Inititate first:
me_iterfig = plt.figure(figsize=(25,8))
me_iterfig_png = png_dir + 'me_iter_sites.png'
me_iterfig_pdf = pdf_dir + 'me_iter_sites.pdf'


# me6cff iterative
plt.errorbar(me6cff_iter.stnum[usta_ind],np.array(me6cff_iter.site_terms)[usta_ind],yerr=np.array(me6cff_iter.site_stderr)[usta_ind],fmt='*',markersize=8,color='#dbb36d',elinewidth=2,capsize=6,capthick=1.5,label='Iterative-Mixed 6 coefficients')

# me5cff iterative
plt.errorbar(me5cff_iter.stnum[usta_ind],np.array(me5cff_iter.site_terms)[usta_ind],yerr=np.array(me5cff_iter.site_stderr)[usta_ind],fmt='p',markersize=8,color='#d17c21',elinewidth=2,capsize=6,capthick=1.5,label='Iterative-Mixed 5 coefficients')


# Station name text:
for stationi in range(len(usta)):
    statext = sm6cff.sta[usta_ind][stationi]
    stx = sm6cff.stnum[usta_ind][stationi]
    sty = -2.1
    rotangle = 90
    
    plt.text(stx,sty,statext,rotation=rotangle,ha='center',va='bottom')
    
plt.plot(me6cff_iter.stnum[usta_ind],(me6cff_iter.site_terms[usta_ind] - me5cff_iter.site_terms[usta_ind]),linestyle='--',color='#7c797a',label='Residual')
plt.plot(me6cff_iter.stnum[usta_ind],me6cff_iter.vs30[usta_ind]/760,linestyle='--',color='#448b8c',label='Vs30/Vref')
plt.plot(me6cff_iter.stnum[usta_ind],(me6cff_iter.site_terms[usta_ind] - me5cff_iter.site_terms[usta_ind])/(me6cff_iter.vs30[usta_ind]/760),linestyle='--',color='#8c445b',label='Residual/Vs30/Vref')

plt.legend(loc=2)

plt.ylim([-2.2,2.5])
plt.title('Site vs. Site term, comparison between iterative ME inversions with 6 and 5 coefficients (Vs30/no Vs30')

plt.savefig(me_iterfig_png)
plt.savefig(me_iterfig_pdf)



## Diff. between path term plots

