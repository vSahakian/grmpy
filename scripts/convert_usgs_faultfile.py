#Make the Anza faults file
#VJS 9/2016

import dread
import cPickle as pickle

#Name of the gmt multisegment file:
multisegpath='/Users/vsahakian/anza/data/faults/Holocene_LatestPleistocene.txt'

#Name of the pckl output file:
#pcklpath='/Users/vsahakian/anza/data/faults/Holocene_LatestPleistocene_117.5w_115.5w_33n_34n.txt'
pcklpath='/Users/vsahakian/anza/data/faults/Holocene_LatestPleistocene_118.0w_115.2w_32.3n_34.5n.pckl'


#Limits to output:
pathlimits=[[-118.0,-115.2],[32.3,34.5]]

#Convert:
allsegments=dread.multiseg2pckl(multisegpath,pcklpath,pathlimits)

#Write to the pickle file:
fout=open(pcklpath,'w')
for segment_i in range(len(allsegments)):
    pickle.dump(allsegments[segment_i],fout)
fout.close()