######Inversion Module######
#VJS 6/2016


import dread as dr
import cdefs as cdf
import numpy as np


#Make the d and G matrices.....

def invinit(db,rng,sdist,smth):
    '''
    Make the d and G matrices for the inversion.  Can use ranges, where the 
    coefficients from teh inversion must be smooth at the edges of the ranges.
    
    Input:
        db:     Database class from cdefs, with:
                mw:     Moment magnitude array (n x 1)
                r:      Distance array (n x 1)
                pga:    log10pga array (n x 1)
        rng:    Array with limits of M ranges, 
                i.e.: [0, 2, 5] means two ranges: M0 - M2, M2 - M5.
        sdist:  Array with distances to include for smoothing (i.e, [1,5,10]
        smth:   Smoothing value
            
    '''
    import numpy as np
    
    print db.mw
    print rng
    

    
    #Number of distances for smoothing:
    smoothnum=len(sdist)
    
    #
    
    
    #Get the indices where each magnitude fits in each bin, digitize index:
    dind=np.digitize(db.mw,rng)
    
    
    
    
        