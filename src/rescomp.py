#Residual Module
#VJS 6/2016


def total_residual(db,d_predicted_ln,vref,predictive_parameter='pga',ncoeff=5,data_correct=-0.6):
    '''
    Compute the total residual and standard deviation for a dataset
    Input:
        db:                      Database object with data to compute
        d_predicted_ln:          Array with predicted value        
        vref:                    Value of reference vs30 
        predictive_parameter:    Predictive parameter. 'pga', or 'pgv'.  Default: 'pga'
        ncoeff:                  Number of coefficients that were inverted for.  Default: 5
        data_correct:            Correct data?  0/correction term = no/by vs30 correction term.  Default: -0.6 correction 
    Output:
        total_residuals:     Array with residual for each data point
        mean_resid:          Mean residual for dataset
        std_dev:             Standard deviation in dataset
    '''
    from numpy import mean,std,log
    
    # Get the data for the parameter being used for the residuals:
    if predictive_parameter=='pga':
        predparam_precorrection=db.pga_pg
    elif predictive_parameter=='pgv':
        predparam_precorrection=db.pgv
        
    #Subtract the prediction from observation... do it in natural log space.
    
    #The PGA from the event object is units of g...not in log10 space.
    if data_correct==0:
        d_observed_ln=log(predparam_precorrection)
        print('\n Using data without any correction \n ')
        print('\n d_observed_ln from predictive parameter is: \n')
        print(d_observed_ln)
    elif ((data_correct!=0) & (predictive_parameter=='pga')):
        vs30correct=data_correct*log(db.vs30/vref)
        d_observed_ln=log(predparam_precorrection)-vs30correct
        print('\n Correcting PGA data by 0.6*ln(vs30/vref) \n')
    else:
        print('You provided PGV data with a vs30 correction...not sure what coefficient to use...not correcting...')
    
    
    print('ln dobserved is ')
    print(d_observed_ln)
    
    print('ln dpredicted is ')
    print(d_predicted_ln)
    
    #However d_predicted is in log10 space...connvert it:
    d_predicted = d_predicted_ln
    d_observed = d_observed_ln
    
    # Get total residual - both of these are already in ln space.
    total_residuals = d_observed - d_predicted   
    
    #Mean total residual?
    mean_resid=mean(total_residuals)
    
    #STandard deviation?
    std_dev=std(total_residuals)
    
    return total_residuals,mean_resid,std_dev
    
def event_residual(eventdb,d_predicted_ln,vref,predictive_parameter='pga',ncoeff=5,data_correct=-0.6):
    '''
    Compute the event residual
    Input:
        eventdb:                 Event object
        d_predicted_ln:          Array with model predictions for each recording in 
                                  event object - predicts ln(PGA)
        vref:                    Reference velocity for vs30 term
        predictive_parameter:    Predictive parameter. 'pga', or 'pgv'.  Default: 'pga'
        ncoeff:                  Number of coefficients that were inverted for.  Default: 5
        data_correct:            Correct data?  0/correction = no/by vs30 correction term.  Default: -0.6 correction 
    Output: 
        E_residual:           Event residual
        E_std_dev:            Standard deviation in the event residual
    '''
    
    from numpy import log,mean,std
    
    #Event number and magnitude?
    evnum=eventdb.evnum[0]
    evmw=eventdb.mw[0]
    
    # Get the data for the parameter being used for the residuals:
    if predictive_parameter=='pga':
        predparam_precorrection=eventdb.pga_pg
    elif predictive_parameter=='pgv':
        predparam_precorrection=eventdb.pgv
    
    #Do all in natural log (np.log) space...
    
    #The PGA from the event object is units of g...not in log10 space.
    #Subtract the predicted value from the event value of pga_pg (in log10 space):
    if data_correct==0:
        d_observed_ln=log(predparam_precorrection)
        print('\n Using data without any correction \n ')
        print('\n d_observed_ln from predictive parameter is: \n')
        print(d_observed_ln)
    elif ((data_correct!=0) & (predictive_parameter=='pga')):
        vs30correct=data_correct*log(eventdb.vs30/vref)
        d_observed_ln=log(predparam_precorrection)-vs30correct
        #print('\n Correcting PGA data by 0.6*ln(vs30/vref) \n'
    else:
        print('You provided PGV data with a vs30 correction...not sure what coefficient to use...exiting')
    
    
    #Get residual for each recording in the event:
    event_residuals=d_observed_ln-d_predicted_ln
    
    #Get the "event residual" (mean of all recordings):
    E_residual=mean(event_residuals)
    E_stderr=std(event_residuals)
    
    #Return these values:
    return evnum,evmw,E_residual, E_stderr
    
def within_event_residual(eventdb,d_predicted_ln,E_residual,vref,predictive_parameter='pga',ncoeff=5,data_correct=-0.6):
    '''
    Compute the within-event residual
    Input:
        eventdb:                 Event object
        d_predicted_ln:          Array with model predictions for each recording in 
                                  event object - predicts ln(PGA)
        E_residual:              Array with event residuals from the event term
        vref:                    Reference velocity for vs30 term
        predictive_parameter:    Predictive parameter. 'pga', or 'pgv'.  Default: 'pga'
        ncoeff:                  Number of coefficients that were inverted for.  Default: 5
        data_correct:            Correct data?  0/correction = no/by vs30 correction term.  Default: -0.6 correction 
        
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
    
    # Get the data for the parameter being used for the residuals:
    if predictive_parameter=='pga':
        predparam_precorrection=eventdb.pga_pg
    elif predictive_parameter=='pgv':
        predparam_precorrection=eventdb.pgv
    
    #Do all in natural log (np.log) space...
    
    #The PGA from the event object is units of g...not in log10 space.
    #Subtract the predicted value from the event value of pga_pg (in log10 space):
    if data_correct==0:
        d_observed_ln=log(predparam_precorrection)
        print('\n Using data without any correction \n ')
        print('\n d_observed_ln from predictive parameter is: \n')
        print(d_observed_ln)
    elif ((data_correct!=0) & (predictive_parameter=='pga')):
        vs30correct=data_correct*log(eventdb.vs30/vref)
        d_observed_ln=log(predparam_precorrection)-vs30correct
        #print('\n Correcting PGA data by -0.6*ln(vs30/vref) \n'
    else:
        print('You provided PGV data with a vs30 correction...not sure what coefficient to use...exiting')
    
    
    #Add the Event residual on to the predicted value, to get the mean 
    #prediction for this event:
    event_predicted=d_predicted_ln+E_residual
    
    #Subtract the value of the prediction for this event from the pga from 
    #each recording in this event
    W_residuals=d_observed_ln-event_predicted
    
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
    