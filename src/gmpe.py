#GMPE
#Module for computing predicted ground motion given coefficents
#VJS 6/2016


def ask2014(db,coeff_file,M1,M2,mdep_ffdf):
    '''
    Compute the predicted ground motionsfor a given set of events using the
    Abrahamson, Silva, and Kamai 2014 model.
    Input:
        db:             Database object with data to plot
        coeff_file:     Path to the file with ASK2014 coefficients
        M1:             M1 referred to in ASK 2014, magnitude scaling break 1
        M2:             M2 referred to in ASK 2014, magnitude scaling break 2
        mdep_ffdf:      Use magnitude dependent fictitous depth?  no=0, yes=1
    Output:
        M:              Magnitude
        PGA:            Predicted PGA
    '''
    
    