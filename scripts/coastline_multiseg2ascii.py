########
## Convert multisegment coastline file from GMT to just points, to plot
##    with scatter in python.
## VJS 4/2018

import numpy as np

## Path:
coastline_multiseg = '/Users/vsahakian/anza/data/faults/socal_coastlines_m.txt'

## Output path:
coastline_ascii = '/Users/vsahakian/anza/data/faults/socal_coastlines.txt'


###########################################3

## Read it in line by line.
f = open(coastline_multiseg,'r')
all_lines = f.readlines()
f.close()


coastline_array_lon = np.array([])
coastline_array_lat = np.array([])

## Parse through them - if it's not started with a carrot, append it:
for iline in all_lines:
    if '>' not in iline:
        i_lon = np.float(iline.split('\n')[0].split('\t')[0])
        i_lat = np.float(iline.split('\n')[0].split('\t')[1])
        
        # Append to coastlien:
        coastline_array_lon = np.r_[coastline_array_lon,i_lon]
        coastline_array_lat = np.r_[coastline_array_lat,i_lat]

coastline_array = np.c_[coastline_array_lon,coastline_array_lat]

## Write to output file:
np.savetxt(coastline_ascii,coastline_array,fmt='%.8f\t%.8f')