#Residual Module
#VJS 6/2016


def total_residual(db,d_predicted_log10):
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
    from numpy import mean,std,log,log10
    
    #Subtract the two... do it in natural log space.
    #The PGA from the event object is %g...not in log10 space.
    pga=db.pga_pg
    #However d_predicted is in log10 space...connvert it:
    d_predicted=10**(d_predicted_log10)
    
    #Now do everything in ln space, since that's what engineers do...
    #ln(pga) - ln(d_predicted).....NOT ln(pga - d_predicted), since
    #in theory d_predicted should predict ln(pga).
    total_residuals=log(pga)-log(d_predicted)
    
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
    