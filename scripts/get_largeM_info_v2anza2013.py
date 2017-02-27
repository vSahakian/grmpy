# Get information about the larger magnitude events that are not
#   recorded at longer distances
#    VJS 2/2017

import cPickle as pickle
import numpy as np

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
    codehome='/home/vsahakian'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    codehome='/Users/vsahakian'

# Residual object path:
dbpath = HOME + '/anza/models/residuals/mixedregr_v2anza2013_Mc_8.5_res4/mixedregr_v2anza2013_Mc_8.5_res4_VR_99.9_robj.pckl'
# Write to an output file:
events_gr3_path = HOME + '/anza/models/residuals/mixedregr_v2anza2013_Mc_8.5_res4/events_mw_greaterthan3.txt'
events_gr35_path = HOME + '/anza/models/residuals/mixedregr_v2anza2013_Mc_8.5_res4/events_mw_greaterthan3_5.txt'

# Open residuals object
rfile = open(dbpath,'r')
robj = pickle.load(rfile)
rfile.close()

# Get magnitude and distance for database
mw = robj.mw
rrup = robj.r

# Find events 
mag_gr3 = np.where(mw > 3)[0]

# Print out:
evnum_gr3 = robj.evnum[mag_gr3]
mw_gr3 = robj.mw[mag_gr3]
ml_gr3 = robj.ml[mag_gr3]
rrup_gr3 = robj.r[mag_gr3]
evlat = robj.elat[mag_gr3]
evlon = robj.elon[mag_gr3]
evdepth = robj.edepth[mag_gr3]
sta = robj.sta[mag_gr3]
stlon = robj.stlon[mag_gr3]
stlat = robj.stlat[mag_gr3]
stelv = robj.stelv[mag_gr3]

# Write to a csv file:
# start writing...
header =  '%5s\t%3s\t%5s\t%5s\t%7s\t%7s\t%5s\t%5s\t%7s\t%7s\t%5s\n' % ('evnum','mw','ml','rrup','elat','elon','edep','sta','stlon','stlat','stelv')
f = open(events_gr3_path,'w')
f.write(header)
for line in range(len(mag_gr3)):
    line_write = '%i\t%.2f\t%.2f\t%.2f\t%.5f\t%.5f\t%.2f\t%5s\t%.5f\t%.5f\t%.2f\n' % (evnum_gr3[line],mw_gr3[line],ml_gr3[line],rrup_gr3[line],evlat[line],evlon[line],evdepth[line],sta[line],stlon[line],stlat[line],stelv[line])
    f.write(line_write)
f.close()

####################################################
####################################################
####################################################

# Open residuals object
rfile = open(dbpath,'r')
robj = pickle.load(rfile)
rfile.close()

# Get magnitude and distance for database
mw = robj.mw
rrup = robj.r

# Find events 
mag_gr = np.where(mw > 3.5)[0]

# Print out:
evnum_gr = robj.evnum[mag_gr]
mw_gr = robj.mw[mag_gr]
ml_gr = robj.ml[mag_gr]
rrup_gr = robj.r[mag_gr]
evlat = robj.elat[mag_gr]
evlon = robj.elon[mag_gr]
evdepth = robj.edepth[mag_gr]
sta = robj.sta[mag_gr]
stlon = robj.stlon[mag_gr]
stlat = robj.stlat[mag_gr]
stelv = robj.stelv[mag_gr]

# Write to a csv file:
# start writing...
header = '%5s\t%3s\t%5s\t%5s\t%7s\t%7s\t%5s\t%5s\t%7s\t%7s\t%5s\n' % ('evnum','mw','ml','rrup','elat','elon','edep','sta','stlon','stlat','stelv')
f = open(events_gr35_path,'w')
f.write(header)
for line in range(len(mag_gr)):
    line_write = '%i\t%.2f\t%.2f\t%.2f\t%.5f\t%.5f\t%.2f\t%5s\t%.5f\t%.5f\t%.2f\n' % (evnum_gr[line],mw_gr[line],ml_gr[line],rrup_gr[line],evlat[line],evlon[line],evdepth[line],sta[line],stlon[line],stlat[line],stelv[line])
   # (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[10])
    f.write(line_write)
f.close()