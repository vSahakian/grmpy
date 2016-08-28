##Test raytracing###
#VJS 8/2016

#Get data info:

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
#run_name='abdb_0-1-2-6.5'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

#home=HOME+'/anza/models/residuals/'
#run_name='abdb_0-6.5_addindex'
#dbpath=HOME+'/anza/data/abdb.pckl'
#modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

home=HOME+'/anza/models/residuals/'
run_name='abdb_5sta_0-6.5_VR'
dbpath=HOME+'/anza/data/abdb_5sta.pckl'
modelpath=HOME+'/anza/models/pckl/regr_0.0_1.0_2.0_6.5_VR_98.9.pckl'

#Get the event and station directories
run_dir=path.expanduser(home+run_name+'/')
eo_dir=run_dir+'event_objs/'
so_dir=run_dir+'sta_objs/'

event_list_path=eo_dir+run_name+'.pckl'
station_list_path=so_dir+run_name+'.pckl'

#Load in db:
dobj=open(dbpath,'r')
db=pickle.load(dobj)
dobj.close()

#Load in list of event and station objects:
eobjs=dread.read_obj_list(event_list_path)
sobjs=dread.read_obj_list(station_list_path)

#Get a test station to use in the residual analysis:
station1=sobjs[0]
station2=sobjs[1]

#Get station info:
sta1=station1.sta
sta2=station2.sta

stnum1=station1.stnum
stnum2=station2.stnum

#Get a few events:
event1num=station1.evnum[0]
event2num=station1.evnum[1]
event3num=station1.evnum[2]

#Looks like event1num is recorded on station1 and station2, as well as event2num.
#check:
station2event1=station2.evnum[0]
station2event3=station2.evnum[1]

if station2event1==event1num:
    keepEvent1=True
    
if station2event3==event3num:
    keepEvent3=True
    
#IF this is true, then events 1, 2, and 3 are recorded on station 1
#   and events 1 and 3 are recorded on station 2.
#This means the source file can contain these three events, and the receiver
#   file can contain stations 1 and 2.  Station 1 will reference all three events;
#       Station 2 will reference event 1 and 3. 

#

#################################
#Ignore the above....
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


