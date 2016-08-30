##Play around iwth the raytracing...
#VJS 8/2016

import run_res
from os import path
import dread
import cPickle as pickle
import raytracing as rt

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=0

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'

#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-6.5_addindex'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

home=HOME+'/anza/models/residuals/'
run_name='abdb_0-6.5_5sta'
dbpath=HOME+'/anza/data/abdb_5sta.pckl'
modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

#Set up the sources.in and receivers.in files for Vs.  So far, Vp has been run:
#in /media/vsahakian/katmai/anza/data/vm/fulltest_Vp

#Run these with raytracing write_sourcein and write_receiverin
#Use the directories defined in lines 15 - 33
#Import raytracing above
#Set velocity type:
veltype=2

#Write source.in file:
rt.write_sourcein(home,run_name,veltype)

#Write recievers.in file:
rt.write_receiverin(home,run_name)


######RUN##########

#For vp:
rayfile=HOME+'/anza/data/vm/fulltest_Vp/rays.dat'
veltype=1

#Read in:
rt.store_rayinfo(home,run_name,rayfile,veltype)

####s

#For vs:
rayfile=HOME+'/anza/data/vm/fulltest_Vs/rays.dat'
veltype=2

#Read in:
rt.store_rayinfo(home,run_name,rayfile,veltype)

###
#Plot rays:
#Vp:
veltype=1
#map:
view=0
axlims=[[-116.9,-116.35],[33.3,33.75]]
stations=1
events=1
by_path=1
mymap='jet'

def plot_rays(home,run_name,veltype,view,axlims,stations,events,by_path,mymap):
    '''
    Plot the raypaths and save the png and pdf figures
    VJS 8/2016
    Input:
        home:           String with the home directory for the project
        run_name:       String with the run name combo of db and model
        veltype:        Vp/Vs (1/2)
        view:           Plot view; map/lat vs depth/lon vs depth, 0/1/2
        axlims:         Axis limits [[xmin, xmax],[ymin,ymax]]
        stations:       Plot stations?  no/yes = 0/1
        events:         Plot events?  no/yes = 0/1
        by_path:        Color raypaths by path?  no/yes=0/1
        mymap:          String with the python colormap to use
        cutoff_val:     Value above/below which to color raypath; otherwise the path is gray (plots if abs(path_term)>=cutoff_val)
    Output: 
        figure          Prints a png and pdf version of the figure to the run fig directory
    '''
    
    from os import path
    
    #Get the run directory:
    run_dir=path.expanduser(home+run_name+'/')
    
    #And the figure directories:
    fig_dir=run_dir+'figs/'
    pdf_dir=fig_dir+'pdfs/'
    
    #Get the residuals object:
    residfile=run_dir+run_name+'_robj.pckl'
    
    #Load the residuals object:
    rfile=open(residfile,'r')
    robj=pickle.load(rfile)
    rfile.close()
    
    figure=robj.plot_raypaths(veltype,view,axlims,stations,events,by_path,mymap)
    
    #Save the figures:
    #Get the figure name:
    pngfile=fig_dir+run_name+'_view'+str(view)+'_sta'+str(stations)+'_ev'+str(events)+'.png'
    pdffile=pdf_dir+run_name+'_view'+str(view)+'_sta'+str(stations)+'_ev'+str(events)+'.png'
    #Save png:
    figure.savefig(pngfile)
    #save pdf:
    figure.savefig(pdffile)