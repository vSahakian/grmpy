######GMPE######
#Module for computing predicted ground motion given coefficents
#VJS 6/2016


def ask2014_pga(M,Rrup,coeff_file,mdep_ffdf,dist_ranges):
    '''
    Compute the predicted ground motionsfor a given set of events using the
    Abrahamson, Silva, and Kamai 2014 model.
    Input:
        M:              Array with moment magnitudes
        Rrup:           Array with Rrup, same size as M
        coeff_file:     Path to the file with ASK2014 coefficients
        mdep_ffdf:      Use magnitude dependent fictitous depth?  no=0, yes=1
        dist_ranges:    
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
    
    
    #Get the PGA coefficients:
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
    def basic(M,Rrup,t_flag):
        '''
        Basic form in computing the predictive paramater (here, PGA) using 
        Abrahamson, Silva, and Kamai 2014's model.
        Input:
            M:          Moment Magnitude
            Rrup:       Closest distance to rupture
            t_flag:     Flag for predictive parameter. 0=PGA.
        Output: 
            f1 (log10(predictive parameter)), 
            to put into full functional form or use alone
        '''
        
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
        #Return the predictive parameter for the basic form, f1:
        return f1,M_sort,f1_sort
    
    
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
    f1,M_sort,f1_sort=basic(M,Rrup,0)
    
    return f1,M_sort,f1_sort