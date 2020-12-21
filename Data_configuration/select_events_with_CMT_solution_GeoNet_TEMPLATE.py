#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:04:50 2020

@author: andrei

This is for locating events in the domain of interest (in space, time, Mw)
"""


from obspy import UTCDateTime
from obspy.clients.fdsn import Client as FDSN_Client
from obspy import read_inventory

client = FDSN_Client("GEONET")

starttime = "2003-11-13 11:00:00.000"
endtime = "2019-11-14 11:00:00.000"
#Get the start and end time for events; lat/lon from VeloModCorners3.txt;
#https://github.com/GeoNet/fdsn/blob/main/examples/GeoNet_FDSN_demo_event.ipynb
cat = client.get_events(starttime=starttime, endtime=endtime,latitude=-42.693,longitude=173.022,maxradius=2.0,minmagnitude=4.8,maxmagnitude=5.5,maxdepth=20)
#cat = client.get_events(starttime=starttime, endtime=endtime,latitude=-42.080445,longitude=173.04739849999999,maxradius=1.4,minmagnitude=4.8,maxmagnitude=5.2,mindepth=10,maxdepth=20)
print(cat)
#print(cat.__str__(print_all=True))
cat.plot(projection="local")

#Create all events' list:
sname = 'events4p8_20km.txt'
fid1=open(sname,'w')
SRF_names = []
for i in range(0, len(cat)):
    fid1.write("%s\n" %(str(cat[i].resource_id).split('/')[1]))
    SRF_names.append(str(cat[i].resource_id).split('/')[1])

#Create all events' list with CMT solution:    
sname_cmt = 'list_'+str(len(SRF_names))+'_4p8_events.csv'
fid2=open(sname_cmt,'w')

#Open all GeoNet CMT database file and search for event in the list "events4p8_20km.txt"    
#https://github.com/GeoNet/data/blob/main/moment-tensor/GeoNet_CMT_solutions.csv
#df = open('GeoNet_CMT_solutions.csv')
df = open('GeoNet_CMT_solutions_20200131.csv')
for i,x in enumerate(df):
    if(i==0):
        #First row in the catalog with EventID, Date,Latitude,Longitude,strike1,dip1,rake1, ...
        fid2.write("%s\n" %(x.strip()))   
    if(i>0):
        row = x.strip().split(',')
        if(row[0] in SRF_names):
            print(row[0])
            fid2.write("%s\n" %(x.strip()))  
            