## Just make Vs30 plots for now...


import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from matplotlib.ticker import MultipleLocator



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


#######################
# Model 3, 6 coefficients, iterative ME with a4, a5, and a6 fixed
me6cff_iter2_path = '/Users/vsahakian/anza/models/residuals/mixedregr_v3anza2013_pga_6coeff_a4_-1.73_a5_-0.01_a6_0.56_Mc_8.5_res4/mixedcoeff_v3anza2013_pga_vs30_6coeff_pga__ncoeff6_Mc_8.5_VR_99.3_a4_-1.73_a5_-0.01_a6_0.56_robj.pckl'

modfile = open(me6cff_iter2_path,'r')
me6cff_iter2 = pickle.load(modfile)
modfile.close()


###########################################
# Unique stations:
usta,usta_ind = np.unique(me6cff_iter.stnum,return_index=True)



###########################################
# Sites vs. vs30....

vs30_5_6fig = plt.figure(figsize=(10,8))
vs30_5_6png = png_dir + 'sites_vs30_5cff_6cff.png'
vs30_5_6pdf = pdf_dir + 'sites_vs30_5cff_6cff.pdf'

# Get correlation values
me6corr,me6r=pearsonr(np.log(me6cff_iter.vs30[usta_ind]),np.array(me6cff_iter.site_terms)[usta_ind])
me5corr,me5r=pearsonr(np.log(me5cff_iter.vs30[usta_ind]),np.array(me5cff_iter.site_terms)[usta_ind])
me6corr2,me6r2=pearsonr(np.log(me6cff_iter2.vs30[usta_ind]),np.array(me6cff_iter2.site_terms)[usta_ind])


# sort by vs30:
vs30sort = np.argsort(me6cff_iter.vs30[usta_ind])

# Get difference to plot:
me_site_diff_unsort = np.array(me6cff_iter2.site_terms)[usta_ind] - np.array(me5cff_iter.site_terms)[usta_ind]
me_site_diff = me_site_diff_unsort[vs30sort]

## Plot..
#plt.scatter(np.log(me6cff_iter.vs30[usta_ind]),np.array(me6cff_iter.site_terms)[usta_ind],s=25,edgecolor='#5d9fba',facecolor='none',marker='D',linewidth=2,label='Mixed iterative, with Vs30, corr coeff = %.2f, p = %.2f' % (me6corr,me6r)) 
#plt.scatter(np.log(me5cff_iter.vs30[usta_ind]),np.array(me5cff_iter.site_terms)[usta_ind],s=25,edgecolor='#b57955',facecolor='none',marker='o',linewidth=2,label='Mixed iterative, without Vs30, corr coeff = %.2f, p = %.2f' % (me5corr, me5r)) 
#plt.scatter(np.log(me6cff_iter2.vs30[usta_ind]),np.array(me6cff_iter2.site_terms)[usta_ind],s=25,edgecolor='#d84936',facecolor='none',marker='^',linewidth=2,label='Mixed iterative, with Vs30 fixed, corr coeff = %.2f, p= %.2f' % (me6corr2,me6r2))


###################
### PLOT 1:  ME6CFF/ME5CFF

# Plot..
plt.errorbar(np.log(me6cff_iter.vs30[usta_ind]),np.array(me6cff_iter.site_terms)[usta_ind],yerr=np.array(me6cff_iter.site_stderr)[usta_ind],markersize=5,color='#5d9fba',fmt='D',elinewidth=2,capsize=4,capthick=1.5,label='Mixed iterative, with Vs30, corr coeff = %.2f, p = %.2f' % (me6corr,me6r)) 
plt.errorbar(np.log(me5cff_iter.vs30[usta_ind]),np.array(me5cff_iter.site_terms)[usta_ind],yerr=np.array(me5cff_iter.site_stderr)[usta_ind],markersize=5,color='#b57955',fmt='o',elinewidth=2,capsize=4,capthick=1.5,label='Mixed iterative, without Vs30, corr coeff = %.2f, p = %.2f' % (me5corr, me5r)) 
#plt.errorbar(np.log(me6cff_iter2.vs30[usta_ind]),np.array(me6cff_iter2.site_terms)[usta_ind],yerr=np.array(me6cff_iter2.site_stderr)[usta_ind],markersize=5,color='#d84936',fmt='^',elinewidth=2,capsize=4,capthick=1.5,label='Mixed iterative, with Vs30 fixed, corr coeff = %.2f, p= %.2f' % (me6corr2,me6r2))

#plt.plot(np.log(me6cff_iter2.vs30[usta_ind][vs30sort]), me_site_diff,linestyle='--',linewidth=2,color='#4f4f4f',label='Residual (vs30 - no vs30)')

ax = vs30_5_6fig.gca()

ax.set_ylim([-2,3.5])

ax.set_xlabel('ln (Vs30) (m/s)')
ax.set_ylabel('Site Term (ln residual)')
ax.set_title('Site term vs. Vs30')

# Set y label spacing:
xmajorLocator = MultipleLocator(1)
ymajorLocator = MultipleLocator(1)

ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_major_locator(xmajorLocator)

# Set y tick spacing
xminorLocator = MultipleLocator(0.1)
yminorLocator = MultipleLocator(0.1)

ax.yaxis.set_minor_locator(yminorLocator)
ax.xaxis.set_minor_locator(xminorLocator)

ax.legend(loc=2)

plt.savefig(vs30_5_6png)
plt.savefig(vs30_5_6pdf)


###################
### PLOT 2:  ME6CFF fixed vs30/ME5CFF

vs30_5_6fix_fig = plt.figure(figsize=(10,8))
vs30_5_6fix_png = png_dir + 'sites_vs30_5cff_6fix.png'
vs30_5_6fix_pdf = pdf_dir + 'sites_vs30_5cff_6fix.pdf'

# Plot..
#plt.errorbar(np.log(me6cff_iter.vs30[usta_ind]),np.array(me6cff_iter.site_terms)[usta_ind],yerr=np.array(me6cff_iter.site_stderr)[usta_ind],markersize=5,color='#5d9fba',fmt='D',elinewidth=2,capsize=4,capthick=1.5,label='Mixed iterative, with Vs30, corr coeff = %.2f, p = %.2f' % (me6corr,me6r)) 
plt.errorbar(np.log(me5cff_iter.vs30[usta_ind]),np.array(me5cff_iter.site_terms)[usta_ind],yerr=np.array(me5cff_iter.site_stderr)[usta_ind],markersize=5,color='#b57955',fmt='o',elinewidth=2,capsize=4,capthick=1.5,label='Mixed iterative, without Vs30, corr coeff = %.2f, p = %.2f' % (me5corr, me5r)) 
plt.errorbar(np.log(me6cff_iter2.vs30[usta_ind]),np.array(me6cff_iter2.site_terms)[usta_ind],yerr=np.array(me6cff_iter2.site_stderr)[usta_ind],markersize=5,color='#d84936',fmt='^',elinewidth=2,capsize=4,capthick=1.5,label='Mixed iterative, with Vs30 fixed, corr coeff = %.2f, p= %.2f' % (me6corr2,me6r2))

#plt.plot(np.log(me6cff_iter2.vs30[usta_ind][vs30sort]), me_site_diff,linestyle='--',linewidth=2,color='#4f4f4f',label='Residual (vs30 - no vs30)')

ax = vs30_5_6fix_fig.gca()

ax.set_ylim([-2,3.5])

ax.set_xlabel('ln (Vs30) (m/s)')
ax.set_ylabel('Site Term (ln residual)')
ax.set_title('Site term vs. Vs30')

# Set y label spacing:
xmajorLocator = MultipleLocator(1)
ymajorLocator = MultipleLocator(1)

ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_major_locator(xmajorLocator)

# Set y tick spacing
xminorLocator = MultipleLocator(0.1)
yminorLocator = MultipleLocator(0.1)

ax.yaxis.set_minor_locator(yminorLocator)
ax.xaxis.set_minor_locator(xminorLocator)

ax.legend(loc=2)

plt.savefig(vs30_5_6fix_png)
plt.savefig(vs30_5_6fix_pdf)


###################
### PLOT 3:  ME6CFF/ME5CFF

vs30_all_fig = plt.figure(figsize=(10,8))
vs30_all_png = png_dir + 'sites_vs30_all.png'
vs30_all_pdf = pdf_dir + 'sites_vs30_all.pdf'

# Plot..
plt.scatter(np.log(me6cff_iter.vs30[usta_ind]),np.array(me6cff_iter.site_terms)[usta_ind],s=25,edgecolor='#5d9fba',facecolor='none',marker='D',linewidth=2,label='Mixed iterative, with Vs30, corr coeff = %.2f, p = %.2f' % (me6corr,me6r)) 
plt.scatter(np.log(me5cff_iter.vs30[usta_ind]),np.array(me5cff_iter.site_terms)[usta_ind],s=25,edgecolor='#b57955',facecolor='none',marker='o',linewidth=2,label='Mixed iterative, without Vs30, corr coeff = %.2f, p = %.2f' % (me5corr, me5r)) 
plt.scatter(np.log(me6cff_iter2.vs30[usta_ind]),np.array(me6cff_iter2.site_terms)[usta_ind],s=25,edgecolor='#d84936',facecolor='none',marker='^',linewidth=2,label='Mixed iterative, with Vs30 fixed, corr coeff = %.2f, p= %.2f' % (me6corr2,me6r2))

ax = vs30_all_fig.gca()

ax.set_ylim([-2,3.5])

ax.set_xlabel('ln (Vs30) (m/s)')
ax.set_ylabel('Site Term (ln residual)')
ax.set_title('Site term vs. Vs30')

# Set y label spacing:
xmajorLocator = MultipleLocator(1)
ymajorLocator = MultipleLocator(1)

ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_major_locator(xmajorLocator)

# Set y tick spacing
xminorLocator = MultipleLocator(0.1)
yminorLocator = MultipleLocator(0.1)

ax.yaxis.set_minor_locator(yminorLocator)
ax.xaxis.set_minor_locator(xminorLocator)

# Set y tick length:
ax.tick_params(which='major',length=7,width=1)
#ax.tick_params(which='minor',length=4,width=1)

ax.legend(loc=2)

plt.savefig(vs30_all_png)
plt.savefig(vs30_all_pdf)