######GMPE######
#Module for computing predicted ground motion given coefficents
#VJS 6/2016


def ask2014_pga(db,coeff_file,M2,mdep_ffdf):
    '''
    Compute the predicted ground motionsfor a given set of events using the
    Abrahamson, Silva, and Kamai 2014 model.
    Input:
        db:             Database object with data to plot
        coeff_file:     Path to the file with ASK2014 coefficients
        M2:             M2 referred to in ASK 2014, magnitude scaling break 2
        mdep_ffdf:      Use magnitude dependent fictitous depth?  no=0, yes=1
    Output:
        M:              Magnitude
        PGA:            Predicted PGA
    '''
    
    from numpy import genfromtxt,where,zeros
    
    #Read in coefficients file:
    ask2014=genfromtxt(coeff_file,skip_header=1)
    
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
    def basic_pga(M,Rrup):
        '''
        Basic form in computing the predictive paramater (here, PGA) using 
        Abrahamson, Silva, and Kamai 2014's model.
        Input:
            M:          Moment Magnitude
            Rrup:       Closest distance to rupture
        Output: 
            f1pga (log10pga), to put into full functional form or use along
        '''
        
        #First, get magnitude dependent fictitious depth, c:
        #Where is it above M5:
        c1_ind=where(M>5)[0]
        c2_ind=where((M<=5) & (M>4))[0]
        c3_ind=where(M<=4)[0]
        
        #Set size of c
        c=zeros(M.size)
        c[c1_ind]=c4[0] - ((c4[0]-1)*(5-M))
        c[c3_ind]=1
        
        #Get geometric spreading distance, R, corrected by c:
        R=(Rrup**2 + c**2)**0.5
        
        ##Compute pga##
        #Depends on the magniutde range, first get the indices for each range:
        m1_ind=where(M>M1[0])[0]
        m2_ind=where((M<M1[0]) & (M>=M2)[0]
        m3_ind=where(M<M2)[0]
        
        #Set the output to zeros, shape of input data (M):
        f1pga=zeros(M.size)
        #Fill in with correct coefficients
        f1pga[m1_ind]=a1 + a5*(M[m1_ind]-M1[0])
    
    
    #Full Functional Form:
    def fmrrup_pga(f1pga,Frv,f7pga,Fn,f8pga,Fas,f11pga,f5pga,Fhw,f4pga,f6pga,f10pga,regional_pga):
        '''
        Compute the predictive parameter (in this case, PGA) using Abrahamson,
        Silva, and Kamai 2014's model. 
        Input: 
        
        Output: 
            log10pga 
        
        '''
        
        
        log10pga = (f1pga + Frv*f7pga + Fn*f8pga + Fas*f11pga + f5pga + \
                        Fhw*f4pga + f6pga + f10pga + regional_pga)
                        
        
        return log10pga
            