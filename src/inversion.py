######Inversion Module######
#VJS 6/2016


import dread as dr
import cdefs as cdf
import numpy as np


#Make the d and G matrices.....

def iinit_pga(db,rng,sdist,smth):
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
    
    #Solving:
    #a1 + a2*M + a3*M^2 + a4*ln(R) + a5*Rrup
    #where R = np.sqrt(R^2 + c^2), 
    #where c=4.5 is "fictitious depth: or "finite fault dimension factor":
    
    #Number of distances for smoothing:
    numsmooth=len(sdist)
    numrng=len(rng)-1
    
    #What are the sizes of G and d:
    #How many data points?
    pgalen=len(db.pga)
    #How many smoothing equations, overall? One set per range...  
    numeq=(numrng*numsmooth)
    #How long will d be then?  Add num of data points and number of smoothing eq
    dlen=pgalen+numeq
    
    #Initiate G and d:
    #d is n x 1, where n was defined above...
    d=np.zeros((dlen,1))
    #G is n x 5*number of ranges:
    G=np.zeros((dlen,numrng*5))
    
    #Get the indices where each magnitude fits in each bin - "digitize index":
    dig_i=np.digitize(db.mw,rng)
    
    #Start filling G and d:
    
    #How many data points in each bin? INitialize...
    numinbin=np.zeros((len(rng)-1,1))
    #Loop over the ranges (len(rng-1) because len(rng) always has one more point
    #than bin)
    for i in range(1,len(rng)):
        #First find how many values in this bin:
        numinbin[i-1]=len(np.where(dig_i==i)[0])
        
    #Loop over each range:
    for j in range(len(rng)-1):
        #Find where the digitize index, dig_i, is equal tot he range we're in:
        bin_i=np.where(dig_i==j+1)[0]
        #Get the magnitudes, distances, and pga's for these indices:
        imw=db.mw[bin_i]
        ir=db.r[bin_i]
        iffdf=db.ffdf[bin_i]
        ipga=db.pga[bin_i]
        
        
    
    

    
    
    
    
        