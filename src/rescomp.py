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
    from numpy import mean,std
    
    #Subtract the two...
    total_residuals=db.pga_pg-d_predicted
    
    #Mean total residual?
    mean_resid=mean(total_residuals)
    
    #STandard deviation?
    std_dev=std(total_residuals)
    
    return total_residuals,mean_resid,std_dev
    
#def event_residual():
#    '''
#    Compute the event residual
#    '''