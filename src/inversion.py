######Inversion Module######
#VJS 6/2016

#Make the d and G matrices.....

def iinit_pga(db,ncoeff,rng,sdist,Mc,smth,mdep_ffdf):
    '''
    Make the d and G matrices for the inversion.  Can use ranges, where the 
    coefficients from teh inversion must be smooth at the edges of the ranges.
    The reason for this is that Baltay and Hanks (2014) showed that different 
    physical processes govern the magnitude-PGA relationship, depending on the 
    magnitude range.  As such, the functional form will not fit the data well 
    for all magnitude ranges, so this allows tailoring of the fit, without 
    adding in the physical processes but while respecting the knowledge that 
    their contributions affect resulting coefficients.
    
    Input:
        db:         Object of class database, from cdefs, with:
                    mw:     Moment magnitude array (n x 1)
                    r:      Distance array (n x 1)
                    pga:    log10pga array (n x 1)
        ncoeff:     Number of coefficients
        rng:        Array with limits of M ranges, 
                    i.e.: [0, 2, 5] means two ranges: M0 - M2, M2 - M5.
        sdist:      Array with distances to include for smoothing (i.e, [1,5,10]
        Mc:         M squared centering term (8.5 in ASK2014)
        smth:       Smoothing value
        mdep_ffdf:  Magnitude-dependent fictitious depth param: 0=off, 1=on
        
            
    '''
    import numpy as np
    
    ######*****************************************************************######
    ###Solving:                                                               ###
    ###                                                                       ###                  
    ###  a1 + a2*M + a3*M^2 + a4*ln(R) + a5*Rrup                              ###
    ###  where R = np.sqrt(R^2 + c^2),                                        ###
    ###  where c is "fictitious depth" or "finite fault dimension factor"     ###
    ###  and is magnitude dependent:                                          ###
    ###        =4.5 for M >5, =4.5 - (4.5 - 1)*(5 - M) for 4<LM<=5,           ###
    ###         =1 for M<=4)                                                  ###
    ######*****************************************************************######
    
    #Get the R for each smoothing distance:
    if mdep_ffdf==0:
        c=4.5
        Rsdist=np.sqrt(sdist**2 + c**2)
        print 'Magnitude dependent fictitious depth is turned OFF'
        
    elif mdep_ffdf==1:
        #Find the indices for each range:
        cr1_ind=np.where(rng>5)[0]
        cr2_ind=np.where((rng<=5) & (rng>4))[0]
        cr3_ind=np.where(rng<=4)[0]
        
        #Zero out the c array, and Rsdist array:
        c=np.zeros(rng.shape)
        Rsdist=np.zeros((len(rng),len(sdist)))
        
        #Find indices of where the magnitudes in rng match the ASK2014 boundaries
        c[cr1_ind]=4.5
        c[cr2_ind]=4.5-((4.5-1)*(5-rng[cr2_ind]))
        c[cr3_ind]=np.ones(cr3_ind.shape)
        
        #Put into Rsdist:
        for i in range(len(sdist)):
            Rsdist[:,i]=np.sqrt(sdist[i]**2 + c**2)
        print 'Magnitude dependent fictitious depth is turned ON'
    
    #Number of distances for smoothing:
    numsmooth=len(sdist)
    numrng=len(rng)-1
    
    #What are the sizes of G and d:
    #How many data points?
    pgalen=len(db.pga_pg)
    
    #How many smoothing equations, overall? One set per range boundary, so do 
    #numrng - 1...  
    numeq=((numrng-1)*numsmooth)
    #How long will d be then?  Add num of data points and number of smoothing eq
    dlen=pgalen+numeq
    
    #Initiate G and d:
    #d is n x 1, where n was defined above...
    d=np.zeros((dlen))
    #G is n x ncoeff*number of ranges:
    G=np.zeros((dlen,(numrng)*ncoeff))
    
    #Get the indices where each magnitude fits in each bin - "digitize index":
    dig_i=np.digitize(db.mw,rng)
    
    
    #####
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
        if mdep_ffdf==0:
            iffdf=db.ffdf[bin_i]
        elif mdep_ffdf==1:
            iffdf=db.md_ffdf[bin_i]
        else:
            print 'Magnitude dependent fictitious depth flag missing.  on = 1, off =0'
        ipga=db.pga_pg[bin_i]
        
        #How many recordings are in this loop/range?
        looplen=numinbin[j]
        
        #With that in mind, populate the G and d matrices with the data, before 
        #the smoothing equations:
        #Name the row and column ranges:
        r_beg=j+crow
        r_end=j+crow+looplen
        c_beg=j+ccol
        c_end=j+ccol+ncoeff
        
        #Parts of G:
        a1=np.ones((numinbin[j]))
        a2=imw
        a3=(Mc-imw)**2
        a4=np.log(iffdf)
        a5=ir
        
        #Define:
        G[r_beg:r_end,c_beg:c_end]=np.c_[a1, a2, a3, a4, a5] 
        d[r_beg:r_end]=np.log10(ipga)
        
        #SMOOTHING:
        #If there are still more ranges after this one, add smoothing so that 
        #the line is continuous between ranges...otherwise, don't add anything.
        if j<(len(rng)-2):  
            print 'Adding smoothing ranges with smoothing %.2f' % (smth)  
            
            #Now fill it with the smoothing info, at the bottom of each range, 
            #except the last range:
            #First, row and column ranges:
            r_beg=j+crow+looplen
            r_end=j+crow+looplen+numsmooth
            c_beg=j+ccol
            c_end=j+ccol+(2*ncoeff)
            
            #Parts of G:
            a1=np.ones((numsmooth))
            a2=np.ones((numsmooth))*rng[j+1]
            a3=np.ones((numsmooth))*(Mc - rng[j+1])**2
            #Add magnitude-dependent fictitious depth:
            if mdep_ffdf==0:
                a4=np.log(Rsdist)
            elif mdep_ffdf==1:
                a4=np.log(Rsdist[j+1])
            a5=sdist
            
            #Define G, d smoothing:
            G[r_beg:r_end, c_beg:c_end]=smth*(np.c_[a1, a2, a3, a4, a5, 
                -1*a1, -1*a2, -1*a3, -1*a4, -1*a5]) 
            d[r_beg:r_end]=np.zeros((numsmooth))
            
        else:
            print 'No smoothing ranges; smoothing not applied'
            
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
    

#def pltGd(G):
#    '''
#    Plot G, d to visualize them. 
#    Input:
#        G:      G matrix for inversion, made by iinit_pga
#        d:      Data matrix for inversion, made by iinit_pga
#    '''
#    
#    import matplotlib.pyplot as plt
    

def invert(G,d):
    '''
    Invert....minimize the L2 norm.
    Input:
        G:      Left hand side matrix
        d:      Data vector
    Output:
        m:      Model vector
    '''
    
    from numpy.linalg import lstsq,norm
    import numpy as np
    
    #Get shape of G and d, check:
    shape_G=G.shape
    shape_d=d.shape
    
    if shape_G[0]!=shape_d[0]:
        if np.size(shape_d)==1:
            print 'G and d dimensions do not agree, G: '+np.str(shape_G[0])+ \
                ' x '+np.str(shape_G[1])+', d: '+np.str(shape_d[0])+' x 1'
        else:
            print 'G and d dimensions do not agree, G: '+np.str(shape_G[0])+ \
                ' x '+np.str(shape_G[1])+', d: '+np.str(shape_d[0])+' x '+ \
                np.str(shape_d[1])
    
    
    #Invert
    m,residual,rank,singular_vals=lstsq(G,d)
    
    #For some reason, lstsq isn't getting residuals for some cases with many 
    #ranges.  #Compute your own:
    
    #L2norm
    L2norm=norm(G.dot(m)-d)
    
    #Get the residual that's supposed to come out of lstsq:
    residual=L2norm**2
    
    #Get the variance reduction:
    VR=(1 - (np.sum(residual)/np.sum(d**2)))*100
    
    print 'L2 norm is %.2f, residual (square) is %.2f, Variance Reduction is \
        %.2f percent' % (L2norm,residual,VR)
    
    
    return m, residual, L2norm, VR, rank, singular_vals
    
        
        
##################################
###Run Mixed Effects Model in R###
##################################

def mixed_effects(pga,m,rrup,vs30,evnum,sta,vref,c,Mc):
    '''
    Run a Mixed effects model to compute the model coefficients (a1 - a5), 
    as well as the event and station terms.  The remaining residual can 
    be classified as the path term plus some aleatory residual.
    
    Input:
        pga:           Array with values of PGA for each recording, in g
        m:             Array with values of moment magnitude per recording
        rrup:          Array with values of Rrup per recording
        vs30:          Array with values of Vs30 per recoridng
        evnum:         Array with values of the event number per recording
        sta:           Array with strings of station names
        vref:          Scalar with value of reference vs30 velocity
        c:             Scalar with fictitious depth parameter (usually 4.5)
        Mc:            Magnitude to center around for M squared functional form component (Mc - M)**2
    Output:
             
    '''
    
    import pandas as pd
    import statsmodels.api as sm
    import numpy as np
    
    ## Set database information
    # Set input for model, that is not "raw" (i.e., M):
    pga_corrected=np.log10(pga) - 0.6*(np.log(vs30/vref))
    m2=(Mc - m)**2
    R=np.sqrt(rrup**2 + c**2)
    lnR=np.log(R)
    
    
    #  First make a dictionary:
    dbdict = {'pga' : pga_corrected, 'm' : m, 'm2' : m2, 'lnR' : lnR, 'rrup' : rrup, 'evnum' : evnum, 'sta' : sta}
    
    # Make datafram ewith Pandas
    data = pd.DataFrame(dbdict)
    
    #Output data to csv:
    data.to_csv('tmp_mixed.csv')

    
    
    #### MAKE SYSTEM CALL TO R ####
    
    
    
    ## Read R results back in
    
    
    
    
    
    
#    ######HISTORY####
#    from rpy2.robjects.package import importr
#from rpy2.robjects.packages import importr
#lme4=importr('lme4')
#import rpy2.robjects as ro
#stats=importr('stats')
#what_home=0
#
#if what_home==0:
#    #Desktop:
#    HOME='/media/vsahakian/katmai'
#elif what_home==1:
#    #Mac:
#    HOME='/Users/vsahakian'
#dbfname=HOME+'/anza/data/databases/db2013_test/db2013test_5sta.pckl'
#dbfname
#import cPickle as pickle
#dbfile=open(dbfname,'r')
#dbfile=open(dbfname,'r')
#db=pickle.load(dbfile)
#dbfile.close()
#M=db.mw
#evnum=db.evnum
#vs30=db.vs30
#sta=db.sta
#Mc=8.1
#pga=db.pga_pg
#rrup=db.r
#
#base=importr('base')
#print(base.R_home())
#print(base._libPaths())
#
#
#
#rpga=ro.FloatVector(pga)
#rm=ro.FloatVector(M)
#
#
#rrup=ro.FloatVector(rrup)
#revnum=ro.FloatVector(evnum)
#
#ro.globalenv["pga"]=pga
#ro.globalenv["pga"]=rpga
#ro.globalenv["m"]=rm
#ro.globalenv["rrup"]=rrup
#ro.globalenv["evnum"]=revnum
#
#testmodel=lme4.lmer("pga ~ m + rrup + (1|evnum)")
#
#
#print(base.summary(testmodel))
#from statsmodels.regression.mixed_linear_model import MixedLMParams
