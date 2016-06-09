######Inversion Module######
#VJS 6/2016

#Make the d and G matrices.....

def iinit_pga(db,rng,sdist,smth):
    '''
    Make the d and G matrices for the inversion.  Can use ranges, where the 
    coefficients from teh inversion must be smooth at the edges of the ranges.
    
    Input:
        db:     Object of class database, from cdefs, with:
                mw:     Moment magnitude array (n x 1)
                r:      Distance array (n x 1)
                pga:    log10pga array (n x 1)
        rng:    Array with limits of M ranges, 
                i.e.: [0, 2, 5] means two ranges: M0 - M2, M2 - M5.
        sdist:  Array with distances to include for smoothing (i.e, [1,5,10]
        smth:   Smoothing value
            
    '''
    import numpy as np
    
    ######*****************************************************************######
    ###Solving:                                                               ###
    ###                                                                       ###                  
    ###  a1 + a2*M + a3*M^2 + a4*ln(R) + a5*Rrup                              ###
    ###  where R = np.sqrt(R^2 + c^2),                                        ###
    ###  where c=4.5 is "fictitious depth" or "finite fault dimension factor" ###
    ######*****************************************************************######
    
    #Get the R for each smoothing distance:
    c=4.5
    Rsdist=np.sqrt(sdist**2 + c**2)
    
    #Number of distances for smoothing:
    numsmooth=len(sdist)
    numrng=len(rng)-1
    
    #What are the sizes of G and d:
    #How many data points?
    pgalen=len(db.pga)
    #How many smoothing equations, overall? One set per range boundary, so do 
    #numrng - 1...  
    numeq=((numrng-1)*numsmooth)
    #How long will d be then?  Add num of data points and number of smoothing eq
    dlen=pgalen+numeq
    
    #Initiate G and d:
    #d is n x 1, where n was defined above...
    d=np.zeros((dlen))
    #G is n x 5*number of ranges:
    G=np.zeros((dlen,(numrng)*5))
    
    #Get the indices where each magnitude fits in each bin - "digitize index":
    dig_i=np.digitize(db.mw,rng)
    
    
    ##Start filling G and d:
    
    #How many data points in each bin? INitialize...
    numinbin=np.zeros((len(rng)-1))
    #Loop over the ranges (len(rng-1) because len(rng) always has one more point
    #than bin)
    for i in range(1,len(rng)):
        #First find how many values in this bin:
        numinbin[i-1]=len(np.where(dig_i==i)[0])
        
    #Convert this to an integer, it is for some reason a float...
    numinbin=numinbin.astype('int')
        
    #Loop over each range:
    #First set a counter for the rows and columns to 0:
    crow=0
    ccol=0
    for j in range(len(rng)-1):
        
        #Find where the digitize index, dig_i, is equal to the range we're in:
        bin_i=np.where(dig_i==j+1)[0]
        
        #Get the magnitudes, distances, and pga's for these indices:
        imw=db.mw[bin_i]
        ir=db.r[bin_i]
        iffdf=db.ffdf[bin_i]
        ipga=db.pga[bin_i]
        
        #How many recordings are in this loop/range?
        looplen=numinbin[j]
        
        #With that in mind, populate the G and d matrices with the data, before 
        #the smoothing equations:
        G[j+crow:j+crow+looplen, j+ccol:j+ccol+5]=np.c_[np.ones((numinbin[j])), 
            imw, (8.5-imw**2), np.log(iffdf), ir] 
        d[(j+crow):(j+crow+looplen)]=np.log10(ipga)
        
        #If there are still more ranges after this one, add smoothing so that 
        #the line is continuous between ranges...otherwise, don't add anything.
        if j<(len(rng)-2):    
            #Now fill it with the smoothing info, at the bottom of each range, 
            #except the last range:
            G[j+crow+looplen:j+crow+looplen+numsmooth, j+ccol:j+ccol+10]=smth*(np.c_[np.ones((numsmooth)), 
                np.ones((numsmooth))*rng[j+1], np.ones((numsmooth))*(8.5 - rng[j+1]**2), 
                np.log(Rsdist), sdist, -1*np.ones((numsmooth)), 
                -1*np.ones((numsmooth))*rng[j+1], -1*np.ones((numsmooth))*(8.5 - rng[j+1]**2), 
                -1*np.log(Rsdist), -1*sdist]) 
            d[j+crow+looplen:j+crow+looplen+numsmooth]=np.zeros((numsmooth))
            
        #To the counter indices, add on:
        #To the rows, add what we are past - so the number of recordings in this
        #range, plus the number of smoothing equations added on, -1 because j
        #increases by 1:
        crow=crow+looplen+numsmooth-1
        #To the columns, add what we are past - we are now once range past, so 
        #we are an extra 5 columns deep (on to the next 5 coefficients), minus
        #one since j increases by 1:
        ccol=ccol+4
    
    
    return G, d
    

    
    
    
    
        