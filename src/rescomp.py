#Residual Module
#VJS 6/2016


def total_residual(db,d_predicted_log10,vref,G,m,d):
    '''
    Compute the total residual and standard deviation for a dataset
    Input:
        db:                  Database object with data to compute
        d_predicted_log10:         Array with predicted value          
    Output:
        total_residuals:     Array with residual for each data point
        mean_resid:          Mean residual for dataset
        std_dev:             Standard deviation in dataset
    '''
    from numpy import mean,std,log,log10,load
    
    #Subtract the two... do it in natural log space.
    #The PGA from the event object is units of g...not in log10 space.
    pga=db.pga_pg
    # Get vs30 and vref:
    vs30correct=0.6*log(db.vs30/vref)
    
    # Get corrected PGA, or d_observed:
    d_observed_log10 = log10(pga) - vs30correct
    
    #However d_predicted is in log10 space...connvert it:
    d_predicted=10**(d_predicted_log10)
    d_observed = 10**(d_observed_log10)
    
    # Predicted residual from inversion:
    inversionmean=mean(G.dot(m)-d)
    print 'Inversion mean is %f' % inversionmean
    
    # Get the difference between d_predicted from G.dot(m) in the inversion and and d_predicted_log10 here (should b the same):
    print 'difference between d_predicted_log10 and G.dot(m) is '
    print (d_predicted_log10 - G.dot(m))
    print '\n difference between d_observed_log10 and d is'
    print (d_observed_log10 - d)
    
    # Load in the inversion data:
    invdat=load('/Users/vsahakian/Desktop/inversiondata.npz')
    # Get those G, m, and d:
    G_inv=invdat['G']
    d_inv=invdat['d']
    m_inv=invdat['m']
    
    #Now get some differences
    print 'Now for loaded inversion info......'
    print 'difference between d_predicted_log10 and inversion G.dot(m) is '
    print (d_predicted_log10 - G_inv.dot(m_inv))
    print '\n difference between d_observed_log10 and inversion d is'
    print (d_observed_log10 - d_inv)
    
    #Now do everything in ln space, since that's what engineers do...
    #ln(pga) - ln(d_predicted).....NOT ln(pga - d_predicted), since
    #in theory d_predicted should predict ln(pga).
    #total_residuals=log(pga)-log(d_predicted)   #-vs30correct
    total_residuals=log(d_observed)-log(d_predicted)   #-vs30correct

    
    #Mean total residual?
    mean_resid=mean(total_residuals)
    
    #STandard deviation?
    std_dev=std(total_residuals)
    
    return total_residuals,mean_resid,std_dev
    
def event_residual(eventdb,d_predicted_log10):
    '''
    Compute the event residual
    Input:
        eventdb:              Event object
        d_predicted_log10:    Array with model predictions for each recording in 
                                  event object
    Output: 
        E_residual:         Event residual
        E_std_dev:          Standard deviation in the event residual
    '''
    
    from numpy import log,log10,mean,std
    
    #Event number and magnitude?
    evnum=eventdb.evnum[0]
    evmw=eventdb.mw[0]
    
    #Do all in natural log (np.log) space...
    
    #The PGA from the event object is %g...not in log10 space.
    #Subtract the predicted value from the event value of pga_pg (in log10 space):
    pga=eventdb.pga_pg
    #Convert d_predicted from its log10 space:
    d_predicted=10**(d_predicted_log10)
    
    #Get residual for each recording in the event:
    event_residuals=log(pga)-log(d_predicted)
    
    #Get the "event residual" (mean of all recordings):
    E_residual=mean(event_residuals)
    E_std_dev=std(event_residuals)
    
    #Return these values:
    return evnum,evmw,E_residual, E_std_dev
    
def within_event_residual(eventdb,d_predicted_log10,E_residual):
    '''
    Compute the within-event residual
    Input:
        
    '''
    
    from numpy import log,log10,mean,std
    
    #Get event number, and magnitude:
    evnum=eventdb.evnum[0]
    evmw=eventdb.mw[0]
    
    #Get the station names and numbers:
    #Name:
    sta=eventdb.sta
    #Number:
    stnum=eventdb.stnum
    
    #Get pga for each event recording - eventdb has pga NOT in log10 space:
    pga=eventdb.pga_pg
    
    #However d_predicted that is piped in IS In log10 space, convert it:
    d_predicted=10**(d_predicted_log10)
    
    #Add the Event residual on to the predicted value, to get the mean 
    #prediction for this event:
    event_predicted=log(d_predicted)+E_residual
    
    #Subtract the value of the prediction for this event from the pga from 
    #each recording in this event
    W_residuals=log(pga)-event_predicted
    
    #Mean and std dev:
    W_mean=mean(W_residuals)
    W_std_dev=std(W_residuals)
    
    return evnum,evmw,sta,stnum,W_residuals,W_mean,W_std_dev
    
    
    
###############################
##### MIXED EFFECTS STUFF #####
###############################

def mixed_event_objects(mixedresiduals):
    '''
    '''
    
    import pandas as pd
    
    # Set mixed_residuals to a shorter name:
    mr=mixedresiduals
    
    # Make a dict structure:
    evdict = {'evnum' : mr.evnum, 'sta' : mr.sta, 'stnum' : mr.stnum, 'ml' : mr.ml, 'mw' : mr.mw, 'pga' : mr.pga, \
        'pgv' : mr.pgv, 'pga_pg' : mr.pga_pg, 'r' : mr.r, 'vs30' : mr.vs30, 'ffdf' : mr.ffdf, 'md_ffdf' : mr.md_ffdf, \
        'elat' : mr.elat, 'elon' : mr.elon, 'edepth' : mr.edepth, 'stlat' : mr.stlat, 'stlon' : mr.stlon, 'stelv' : mr.stelv, \
        'source_i' : mr.source_i, 'receiver_i' : mr.receiver_i}     

    allevents = pd.DataFrame(evdict)
    
    allevents[allevents.evnum == 668267]
    