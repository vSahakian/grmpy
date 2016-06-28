#Residual Module
#VJS 6/2016


def total_residual(db,d_predicted):
    '''
    Compute the total residual and standard deviation for a dataset
    Input:
        db:                  Database object with data to compute
        d_predicted:         Array with predicted value          
    Output:
        total_residuals:     Array with residual for each data point
        mean_resid:          Mean residual for dataset
        std_dev:             Standard deviation in dataset
    '''
    from numpy import mean,std,log10
    
    #Subtract the two...
    pga=log10(db.pga_pg)
    total_residuals=pga-d_predicted
    
    #Mean total residual?
    mean_resid=mean(total_residuals)
    
    #STandard deviation?
    std_dev=std(total_residuals)
    
    return total_residuals,mean_resid,std_dev
    
def event_residual(eventdb,d_predicted):
    '''
    Compute the event residual
    Input:
        eventdb:            Event object
        d_predicted:        Array with model predictions for each recording in 
                            event object
    Output: 
        E_residual:         Event residual
        E_std_dev:          Standard deviation in the event residual
    '''
    
    from numpy import log10,mean,std
    
    #Subtract the predicted value from the event value of pga_pg (in log10 space):
    pga=log10(eventdb.pga_pg)
    #Get residual for each recording in the event:
    event_residuals=pga-d_predicted
    
    #Get the "event residual" (mean of all recordings):
    E_residual=mean(event_residuals)
    E_std_dev=std(event_residuals)
    
    #Return these values:
    return E_residual, E_std_dev
    
def within_event_residual(eventdb,d_predicted,E_residual):
    '''
    Compute the within-event residual
    Input:
        
    '''
    
    from numpy import log10,mean,std
    
    #Get pga for each event recording:
    pga=log10(eventdb.pga_pg)
    
    #Add the Event residual on to the predicted value, to get the mean 
    #prediction for this event:
    event_predicted=d_predicted+E_residual
    
    #Subtract the value of the prediction for this event from the pga from 
    #each recording in this event
    W_residuals=pga-event_predicted
    
    #Mean and std dev:
    
    
    