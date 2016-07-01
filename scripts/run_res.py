####Run the residual analysis####
'''
Run the Residual analysis for a given parameter file
VJS 7/2016
'''


#Initialize the folders
def init(home,run_name):
    '''
    Initialize the folders for output for a new run (new database/model combo)
    Input:
        home:           String to path of the home environment
        run_name:       String with the name for this db/model combo
    Output: 
        Makes directories...
    '''
    
    from shutil import rmtree,copy
    from os import path,makedirs,environ
    
    clob='y'
    run_dir=path.expanduser(home+run_name+'/')
    if path.exists(run_dir):    #If the path exists...
        print 'Run directory is: '+run_dir
        clob=raw_input('The run directory already exists, are you sure you want to clobber it? (y/n)')
        if clob is 'y' or clob is 'Y':
            clob=raw_input('This deletes everything in '+run_dir+', are you sure you want to do that? (y/n)')
            if clob is 'y' or clob is 'Y':
                rmtree(run_dir)
            else:
                #Don't do anything
                print 'Aborting mission...not clobbering...'
        else:
            #AGain, don't do anything
            print 'Aborting mission...not clobbering...'
    if clob is 'y' or clob is 'Y':
        #Make the main directory
        makedirs(run_dir)
        #And the subdirectories:
        makedirs(run_dir+'event_objs/')
        makedirs(run_dir+'sta_objs/')
        makedirs(run_dir+'E_resids/')
        makedirs(run_dir+'W_resids/')
        makedirs(run_dir+'site_resids/')
        makedirs(run_dir+'figs/')
        
        