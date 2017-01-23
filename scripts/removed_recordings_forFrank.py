##
# Read in the pickle file with removed recordings, and get the original database information from it
##
# VJS 1/2017

import cPickle as pickle


# Residuals object path:
robjpath = '/media/vsahakian/katmai/anza/models/residuals/v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj_rejected.pckl'

#######################################################
#######################################################
#######################################################
# This reads in information for every recording that was removed (so some events/stations may be repeated)

# Opening the object:
removed_file=open(robjpath,'r')
removed_obj=pickle.load(removed_file)
removed_file.close()

# Event/source info:
evnum = removed_obj.evnum
evlat = removed_obj.elat
evlon = removed_obj.elon
evdepth = removed_obj.edepth
mw = removed_obj.mw
ml = removed_obj.ml

# Station info:
sta = removed_obj.sta
stnum = removed_obj.stnum
stlat = removed_obj.stlat
stlon = removed_obj.stlon
stelv = removed_obj.stelv

# Other info:
pga = removed_obj.pga     # PGA
pgv = removed_obj.pgv     # PGV
rrup = removed_obj.r
total_residual = removed_obj.total_residual     # Total residual



