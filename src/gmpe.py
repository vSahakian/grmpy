######GMPE######
#Module for computing predicted ground motion given coefficents
#VJS 6/2016


def compute_model(m,rng,mw,r,ffdf,vs30,Mc,vref,mdep_ffdf,ncoeff):
    '''
    Compute the predicted value for a certain mw and Rrup 
    Input:
        m:              Model parameters/coefficients resulting from inversion
        rng:            Magnitude ranges used in inversion
        mw:             Array of moment magnitudes for which to compute
        r:              Array of distances associated with the mw
        ffdf:           Array with the fictitious depth R - see flag
        vs30:           Array with vs30 values
        Mc:             Value to center squared magnitude term around (i.e., 8.5 or 8.1, etc.)
        vref:           Vref value
        mdep_ffdf:      Magnitude dependent fictitious depth flag
        ncoeff:         Number of coefficients inverted for.
    Output:
        d_predicted:    Predicted value of ln(predictive parameter) at that data point    
    '''
    import numpy as np
    
    #Find which magnitude range each mw belongs in, to know which set of 
    #coefficients to use to compute it:
    dig_i=np.digitize(mw,rng)
    
    #Zero out the output d vector, and append to it later on:
    d_predicted=np.array([])
    
    #For each range, get the coefficients:
    for range_i in range(len(rng)-1):
        #Get the coefficients of the model for this magnitude range:
        a1=m[range_i*ncoeff]
        a2=m[(range_i*ncoeff)+1]
        a3=m[(range_i*ncoeff)+2]
        a4=m[(range_i*ncoeff)+3]
        a5=m[(range_i*ncoeff)+4]
        
        if ncoeff==6:
            a6=m[(range_i*ncoeff)+5]
    
        #Where is mw in this range?
        bin_i=np.where(dig_i==range_i+1)
        
        #What input data is needed to compute?
        #Mag, R (ffdf), and vs30...
        mw_rangei=mw[bin_i]
        r_rangei=r[bin_i]
        ffdf_rangei=ffdf[bin_i]
        vs30_rangei=vs30[bin_i]
        
        #Now compute the predicted value:    
        if ncoeff==5:
            d_predicted_i=a1+a2*mw_rangei + a3*(Mc-mw_rangei)**2 + a4*np.log(ffdf_rangei) + \
                a5*r_rangei 
            print 'predicting with just 5 coefficients...'
            print '%f, %f, %f, %f, %f ' % (a1, a2, a3, a4, a5)
        elif ncoeff==6:
            d_predicted_i=a1+a2*mw_rangei + a3*(Mc-mw_rangei)**2 + a4*np.log(ffdf_rangei) + \
                a5*r_rangei + a6*np.log(vs30_rangei/vref) 
            print 'predicting with 6 coefficients...'
            print '%f, %f, %f, %f, %f, %f ' % (a1, a2, a3, a4, a5, a6)
        
        #d_predicted_i=a1+a2*mw_rangei + a3*(Mc-mw_rangei)**2 + a4*np.log(ffdf_rangei) + \
        #        a5*r_rangei + 0.6*np.log(vs30_rangei/vref)  
            
        ######^^NOT USING ANYMMORE.....CREATES A BIAS IN RESIDUALS!!  Instead, remove before inverting... so ln(pga) - 0.6*ln(vs30/vref)###
        
        #And append to the final output value:
        d_predicted=np.r_[d_predicted,d_predicted_i]
        
    return d_predicted



###############################################################################
def compute_model_fixeddist(m,rng,sdist,Mc,mdep_ffdf,ncoeff=5):
    '''
    Compute the model values given the coefficients resulting from an inversion;
    Obtain values for given distances.
    
    Input:
        m:          Model coefficients
        rng:        Array with ranges used to get the model coefficients
        sdist:      Array with distances at which to compute the model
        Mc:         Magnitude squred term (i.e., 8.5 or 8.1)
        mdep_ffdf:  Flag for fictitious depth mag-dependence; 0=no, 1=yes
        ncoeff:     Number of coefficients inverted for. Default: 5
    Output:
        mw_out:     Array of mw
        d_out:      Predicted ln(PGA)
    '''
    
    import numpy as np
    
    ####
    #Magnitude dependence?
    if mdep_ffdf==0:
        print 'Magnitude dependent fictitious depth is OFF - check you provided the right ffdf'
    elif mdep_ffdf==1:
        print 'Magnitude dependent fictitious depth is ON - check you provided the right ffdf'
    else:
        print 'Magnitude dependent fictitous depth flag not provided correctly: OFF=0, ON=1'
            
            
    #Loop over distances and ranges to get model output for each distance range, R.
    
    #Loop over distances, get the values for each distance first:
    for j in range(len(sdist)):
        #First, make empty arrays to append computed model/magnitude to:
        mw_range=np.array([])
        d_range=np.array([])
        
        #Then Loop over the ranges, get data for each range:
        for i in range(len(rng)-1):                                   
            #Get the magnitudes to plot against:
            mw=np.linspace(rng[i],rng[i+1],100)
            
            #Get magnitude dependent fictitious depth??
            #Fictitous depth coefficient:
            c4=4.5
            
            #If it's not agnitude dependent, use the scalar coefficient; if not,
            #Use the rules in ASK 2014
            if mdep_ffdf==0:
                R=np.sqrt(sdist[j]**2 + c4**2)
            elif mdep_ffdf==1:
                #Find the indices for each range:
                cr1_ind=np.where(mw>5)
                cr2_ind=np.where((mw<=5) & (mw>4))
                cr3_ind=np.where(mw<=4)
                
                #Zero out the c array:
                c=np.zeros(mw.shape)
                c[cr1_ind]=c4
                c[cr2_ind]=c4-((c4-1)*(5-mw[cr2_ind]))
                c[cr3_ind]=1
                R=np.sqrt(sdist[j]**2 + c**2)
            
            #Get the coefficients for this range:
            a1=m[i*ncoeff]
            a2=m[(i*ncoeff)+1]
            a3=m[(i*ncoeff)+2]
            a4=m[(i*ncoeff)+3]
            a5=m[(i*ncoeff)+4]
        
            
            #    GMPE:    #
            # This now predicts lnPGA, but it goes into somethign that plots log10PGA...
            #       So convert to log10...
            d_ln=a1+a2*mw + a3*(Mc-mw)**2 + a4*np.log(R) + \
                a5*sdist[j] 

            
            #Add these onto the bigger array, for this range:     
            mw_range=np.r_[mw_range,mw]
            d_range=np.r_[d_range,d_ln]
        
        #Now add these onto the final arrayt, horizontally:
        if j==0:
            mw_out=mw_range
            d_out=d_range
        else:
            mw_out=np.c_[mw_out,mw_range]
            d_out=np.c_[d_out,d_range]

    #Print them out...
    return mw_out,d_out


###############################################################################
def ask2014(M,Rrup,coeff_file,mdep_ffdf,dist_ranges,ncoeff=5,predictive_parameter='pga'):
    '''
    Compute the predicted ground motionsfor a given set of events using the
    Abrahamson, Silva, and Kamai 2014 model.
    Input:
        M:                      Array with moment magnitudes
        Rrup:                   Array with Rrup, same size as M
        coeff_file:             Path to the file with ASK2014 coefficients
        mdep_ffdf:              Use magnitude dependent fictitous depth?  no=0, yes=1
        dist_ranges:    
        ncoeff:                 Number of coefficients in inversion.  Default: 5
        predictive_parameter:   Predictive parameter: 'pga', or 'pgv'.  Default: 'pga'
    Output:
        M:              Magnitude
        PGA:            Predicted PGA
    '''
    
    from numpy import genfromtxt,where,zeros,log,argsort
    
    #Read in coefficients file:
    ask2014=genfromtxt(coeff_file,skip_header=1)
    
    #Define period-independent constants:
    M2=5   #Magnitude scaling break #2 - #1 is period dependent
    a7=0
    n=1.5
    
    
    #Get the coefficients:
    #Coefficients for median ground motions
    T=ask2014[:,0]
    M1=ask2014[:,1]
    vlin=ask2014[:,2]
    b=ask2014[:,3]
    c=ask2014[:,4]
    c4=ask2014[:,5]
    a1=ask2014[:,6]
    a2=ask2014[:,7]
    a3=ask2014[:,8]
    a4=ask2014[:,9]
    a5=ask2014[:,10]
    a6=ask2014[:,11]
    a8=ask2014[:,12]
    a10=ask2014[:,13]
    a11=ask2014[:,14]
    a12=ask2014[:,15]
    a13=ask2014[:,16]
    a14=ask2014[:,17]
    a15=ask2014[:,18]
    a17=ask2014[:,19]
    a43=ask2014[:,20]
    a44=ask2014[:,21]
    a45=ask2014[:,22]
    a46=ask2014[:,23]
    
    #Coefficients for the median ground motion in other regions
    a25=ask2014[:,24]
    a28=ask2014[:,25]
    a29=ask2014[:,26]
    a31=ask2014[:,27]
    a36=ask2014[:,28]
    a37=ask2014[:,29]
    a38=ask2014[:,30]
    a39=ask2014[:,31]
    a40=ask2014[:,32]
    a41=ask2014[:,33]
    a42=ask2014[:,34]
    
    #Coefficients for hte standard deviation
    s1e=ask2014[:,35]
    s2e=ask2014[:,36]
    s3=ask2014[:,37]
    s4=ask2014[:,38]
    s1m=ask2014[:,39]
    s2m=ask2014[:,40]
    s5=ask2014[:,41]
    s6=ask2014[:,42]
    
    
    ### Model ###
    #Full Model:
    #f(M,Rrup) = (f1pga + Frv*f7pga + Fn*f8pga + Fas*f11pga + f5pga + Fhw*f4pga + f6pga + f10pga + regional_pga)
    #Frv, Fn, Fas, Fhw are flags to turn on/off reverse faulting, normal faulting, and aftershocks, respectively.
    
    #Basic model:
    def basic(M,Rrup,t_flag,ncoeff):
        '''
        Basic form in computing the predictive paramater (here, PGA) using 
        Abrahamson, Silva, and Kamai 2014's model.
        Input:
            M:          Moment Magnitude
            Rrup:       Closest distance to rupture
            t_flag:     Flag for predictive parameter. 0=PGA, 1=PGV
            ncoeff:     Number of coefficients in version
        Output: 
            f1 (log10(predictive parameter)), 
            to put into full functional form or use alone
        '''
        
        from numpy import log10,exp
        
        #Define coefficients for given predictive parameter:
        M1t=M1[t_flag]
        c4t=c4[t_flag]
        a1t=a1[t_flag]
        a2t=a2[t_flag]
        a3t=a3[t_flag]
        a4t=a4[t_flag]
        a5t=a5[t_flag]
        a6t=a6[t_flag]   
        a8t=a8[t_flag]
        a17t=a17[t_flag]
        
        ##Period-independent:
        a7t=a7 
        M2t=M2
        
        
        #First, get magnitude dependent fictitious depth, c:
        #Where is it above M5:
        c1_ind=where(M>5)[0]
        c2_ind=where((M<=5) & (M>4))[0]
        c3_ind=where(M<=4)[0]
        
        #Set size of c
        c=zeros(M.size)
        c[c1_ind]=c4t
        c[c2_ind]=c4t - ((c4t-1)*(5-M[c2_ind]))
        c[c3_ind]=1
        
        #Get geometric spreading distance, R, corrected by c:
        R=(Rrup**2 + c**2)**0.5
        
        ##Compute pga##
        #Depends on the magniutde range, first get the indices for each range:
        m1_ind=where(M>M1t)[0]
        m2_ind=where((M<M1t) & (M>=M2))[0]
        m3_ind=where(M<M2)[0]
        
        #Set the output to zeros, shape of input data (M):
        f1=zeros(M.size)
        #Fill in with correct coefficients
        f1[m1_ind] = a1t + a5t*(M[m1_ind]-M1t) + \
                        a8t*(8.5 - M[m1_ind])**2 + \
                        (a2t + a3t*(M[m1_ind] - M1t))*log(R[m1_ind]) + \
                        a17t*Rrup[m1_ind]
        
        f1[m2_ind] = a1t + a4t*(M[m2_ind] - M1t) + a8t*(8.5 - M[m2_ind])**2 + \
                        (a2t + a3t*(M[m2_ind] - M1t))*log(R[m2_ind]) + \
                        a17t*Rrup[m2_ind]
                        
        f1[m3_ind] = a1t + a4t*(M2t - M1t) + a8t*(8.5 - M2)**2 + \
                        a6t*(M[m3_ind] - M2t) + a7t*(M[m3_ind] -M2t)**2 + \
                        (a2t + a3t*(M2t - M1t))*log(R[m3_ind]) + \
                        a17t*Rrup[m3_ind]
                        
                        
        #Sort them by magnitude also, for plotting:
        sort_ind=argsort(M)
        M_sort=M[sort_ind]
        f1_sort=f1[sort_ind]
        
        #Convert them to log10 space....right now they're in ln space...
        f1_sort_log10=log10(exp(f1_sort))
        
        #Return the predictive parameter for the basic form, f1:
        return f1,M_sort,f1_sort_log10
    
    
    # #Full Functional Form:
    #def fmrrup(f1,Frv,f7,Fn,f8,Fas,f11,f5,Fhw,f4,f6,f10,regional,t_flag):
    #    '''
    #    Compute the predictive parameter (in this case, PGA) using Abrahamson,
    #    Silva, and Kamai 2014's model. 
    #    Input: 
    #    
    #    Output: 
    #        log10pga 
    #    
    #    '''
    #    
    #    
    #    #log10pga = (f1 + Frv*f7 + Fn*f8 + Fas*f11 + f5 + \
    #                    Fhw*f4 + f6 + f10 + regional)
    #                    
    #    
    #    #return log10pp
    #        
            
    
    ####TEMPORARY....JUST FOR GETTING THE BASIC FORM...###
    if predictive_parameter=='pga':
        #Set t_flag to 0 (pga)
        f1,M_sort,f1_sort=basic(M,Rrup,0,ncoeff)
    elif predictive_parameter=='pgv':
        f1,M_sort,f1_sort=basic(M,Rrup,1,ncoeff)
    
    return f1,M_sort,f1_sort