## Add to Janine's data file a column that has a flag for whether or not
##      that recording was kept or rejected

# Read in the residuals object
robjpath = '/Users/vsahakian/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj_rejected.pckl'

# Also the flatfile from Janine:
jsb_flatfile = '/Users/vsahakian/anza/data/databases/v2anza2013/PGA_2013_MW_allC_qtime.dat'

# Output file to write to for Janine:
jsb_rejrecordings_flatfile = '/Users/vsahakian/anza/data/databases/v2anza2013/PGA_2013_MW_allC_qtime_removalflag.dat'

################################################

from numpy import genfromtxt,zeros,where
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
dat_r=genfromtxt(jsb_flatfile,usecols=range(2,15))

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
origin_time_r=dat_r[:,12]


####
## Set up new structure to write to file..  ##

# Zero out an array, to use for a 'rejection' flag
rejection_flag=zeros(len(mw_r))

# Loop over the recordings from the flatfile.  If that recording matches one of hte rejected recordings, 
for irecord in range(len(rejection_flag)):
    evnum_recordi = evnum_r[irecord]
    sta_recordi = sta_r[irecord]
    
    # Where in the database of rejected recordings is this record/does it exist?
    rejected_index = where((robj.evnum==evnum_recordi) & (robj.sta==sta_recordi))[0]
    
    # If it exists/was rejected, set rejection_flag==1:
    if rejected_index:
        rejection_flag[irecord] = 1


# Now print them all out to a file:
# Header:
header = '%5s\t%3s\t%9s\t%11s\t%8s\t%6s\t%9s\t%11s\t%9s\t%8s\t%11s\t%8s\t%8s\t%10s\t%17s\t%10s\t%s\n' % ('sta', 'chan', 'stlat', 'stlon', 'stelv', 'evnum', 'evlat', 'evlon', 'evdepth', 'grcircle', 'ml', 'mw', 'pga', 'pga_snr', 'origintime', 'date', 'removal_flag')
#open file to write:
f = open(jsb_rejrecordings_flatfile,'w')
f.write(header)
# Loop over lines and write:
for linei in range(len(mw_r)):
    linewrite = '%5s\t%3s\t%9.6f\t%11.6f\t%8.6f\t%6i\t%9.6f\t%11.6f\t%9.6f\t%8.6f\t%11.6f\t%8.6f\t%8.6f\t%10.6f\t%17.6f\t%10s\t%i\n' % (sta_r[linei], chan_r[linei], stlat_r[linei], stlon_r[linei], stelv_r[linei], evnum_r[linei], evlat_r[linei], evlon_r[linei], evdep_r[linei], grcircle_r[linei], ml_r[linei], mw_r[linei], pga_mgal_r[linei], pga_sn_r[linei], origin_time_r[linei], origindate_r[linei], rejection_flag[linei])
    f.write(linewrite)
f.close()
