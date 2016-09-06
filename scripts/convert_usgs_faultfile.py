#Make the Anza faults file
#VJS 9/2016

import dread

#Name of the gmt multisegment file:
multisegpath='/Users/vsahakian/anza/data/faults/Holocene_LatestPleistocene.txt'

#Name of the pckl output file:
pcklpath='/Users/vsahakian/anza/data/faults/Holocene_LatestPleistocene_117.5w_115.5w_33n_34n.txt'

#Limits to output:
pathlimits=[[-117.5,-115.5],[33,34]]

#Convert:
allsegments=dread.multiseg2pckl(multisegpath,pcklpath,pathlimits)