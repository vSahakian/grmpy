## Add to Janine's data file a column that has a flag for whether or not
##      that recording was kept or rejected

# Read in the residuals object
robjpath = '/Users/vsahakian/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj_rejected.pckl'

# Also the flatfile from Janine:
jsb_flatfile = '/Users/vsahakian/anza/data/databases/v2anza2013/PGA_2013_MW_allC_qtime.dat'

################################################

from numpy import genfromtxt,zeros
import cPickle as pickle

## Read in rejected residuals object:
rfile=open(robjpath,'r')
robj=pickle.load(rfile)
rfile.close()

## Read in jsb file:
# First stations and channel
sta_r=genfromtxt(jsb_flatfile,dtype='S5',usecols=[0])
chan_r=genfromtxt(jsb_flatfile,dtype='S5',usecols=[1])
origindate_r=genfromtxt(jsb_flatfile,dtype='S10',usecols=[15], delimiter=' ')

#Then the float data from the flatfile:
dat_r=genfromtxt(jsb_flatfile,usecols=range(2,14))

#Set variables from dat_r:
stlat_r=dat_r[:,0]
stlon_r=dat_r[:,1]
stelv_r=dat_r[:,2]
evnum_r=dat_r[:,3]
evlat_r=dat_r[:,4]
evlon_r=dat_r[:,5]
evdep_r=dat_r[:,6]
grcircle_r=dat_r[:,7]
ml_r=dat_r[:,8]
mw_r=dat_r[:,9]
pga_mgal_r=dat_r[:,10]
pga_sn_r=dat_r[:,11]
origin_time_r=dat_r[:,13]


####
## Set up new structure to write to file..  ##

# Zero out an array, to use for a 'rejection' flag
rejection_flag=zeros(len(mw_r))

# Loop over the recordings from the flatfile.  If that recording matches one of hte rejected recordings, 

