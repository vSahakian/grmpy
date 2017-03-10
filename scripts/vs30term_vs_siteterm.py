# Make a plot of the site terms vs vs30 terms to look at what's going on with
#     vs30 coefficient in the inversion

import cPickle as pickle
from numpy import unique,array
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# database path:
dbpath='/Users/vsahakian/anza/models/residuals/v2anza2013_pga_no_vs30_Mc8.5_pgrid_5sta_res4/v2anza2013_pga_no_vs30_Mc8.5_pgrid_5sta_res4_robj.pckl'
figfile='/Users/vsahakian/Desktop/vs30_vs_siteterm.png'


# Read in database:
dbfile=open(dbpath,'r')
db=pickle.load(dbfile)
dbfile.close()

# Indices of unique stations:
unique_sta,unique_sta_ind=unique(db.sta,return_index=True)

# Unique station site terms:
unique_vs30 = db.vs30[unique_sta_ind]
unique_siteterm = array(db.site_terms)[unique_sta_ind]

#Pearsons coefficient
pcoeff = pearsonr(unique_vs30,unique_siteterm)

plt.scatter(unique_vs30,unique_siteterm,edgecolors='k',facecolors='none',linewidth=2)
plt.xlabel('Vs30 (m/s)')
plt.ylabel('Site term (ln residual)')
plt.title('Vs30 vs Single-mean Site term for Unique stations, ccoeff = %.2f, p = %.2f' % (pcoeff[0],pcoeff[1]))
plt.show()
plt.savefig(figfile)