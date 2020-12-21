#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 12:48:10 2020

@author: andrei
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def symmetric_source(SRF_all,m_lon,m_lat,a,b):
    nShot = np.shape(SRF_all)[1]
    SRF_index = np.array(SRF_all[0,:])
#    nShot_UL = int(np.sum(SRF_all[0,:]))
    for i in range(0,nShot):
        if(SRF_all[0,i] == 1):
            #symmetric point cross center
            sym_lon = 2*m_lon - SRF_all[1,i]                
            sym_lat = 2*m_lat - SRF_all[2,i]
            sym_distance = np.square(sym_lon*np.ones(nShot)-SRF_all[1,:])+np.square(sym_lat*np.ones(nShot)-SRF_all[2,:])
            sym_max = np.max(sym_distance)
            for j in range(0,nShot):
                if(SRF_all[0,j] == 1):
                    sym_distance[j] = sym_max
            sym_id = np.argmin(sym_distance)            
            SRF_index[sym_id] = 1

            #symmetric point cross diagonal            
            i_lon = (SRF_all[1,i]+(SRF_all[2,i]-b)/a)/2
            i_lat = (SRF_all[2,i]+(a*SRF_all[1,i]+b))/2
            
            sym_lon = 2*i_lon - SRF_all[1,i]                
            sym_lat = 2*i_lat - SRF_all[2,i]
            sym_distance = np.square(sym_lon*np.ones(nShot)-SRF_all[1,:])+np.square(sym_lat*np.ones(nShot)-SRF_all[2,:])
            sym_max = np.max(sym_distance)
            for j in range(0,nShot):
                if(SRF_all[0,j] == 1):
                    sym_distance[j] = sym_max
            sym_id = np.argmin(sym_distance)            
            SRF_index[sym_id] = 1            
            
    return SRF_index
                
def read_stat_name_utm(station_file):

    with open(station_file, 'r') as f:
        lines = f.readlines()
    line0=lines[0].split()
    nRec=int(line0[0])
    R=np.zeros((nRec,2))
    statnames = [] 
    for i in range(1,nRec+1):
        line_i=lines[i].split()
        R[i-1,0]=float(line_i[3])
        R[i-1,1]=float(line_i[4])
        statnames.append(line_i[2])
    return nRec, R, statnames                
############################
plt.figure(figsize=(5,5))

station_file = 'STATION_utm_dh_4km.txt' 
nRec, R, statnames = read_stat_name_utm(station_file)    

nlon,nlat = np.loadtxt('North_WGS84.txt',unpack=True)
slon,slat = np.loadtxt('South_WGS84.txt',unpack=True)

#Get lon/lat manually
lonmin = 171
lonmax = 175
latmin = -43
latmax = -40

lon = [lonmin,lonmax,lonmax,lonmin,lonmin]
lat = [latmax,latmax,latmin,latmin,latmax]

#Get lon/lat from model file
lon = []
lat = []
f = open('VeloModCorners3.txt')
for i,x in enumerate(f):
    if i > 1:
        row = x.strip().split()
        lon.append(float(row[0]))
        lat.append(float(row[1]))
lon.append(lon[0])
lat.append(lat[0])

plt.plot(nlon,nlat,c='grey')
plt.plot(slon,slat,c='grey')

plt.plot([lon[0], lon[2]],[lat[0], lat[2]],c='r')
##y = ax+b
#a = (lat[0]-lat[2])/(lon[0]-lon[2])
#b = lat[0] - a*lon[0]

lat_u = (lat[2]+lat[1])/2
lon_u = (lon[2]+lon[1])/2

lat_d = (lat[0]+lat[3])/2
lon_d = (lon[0]+lon[3])/2

#plt.plot([lon_u, lon_d],[lat_u, lat_d],c='r')
a = (lat_u-lat_d)/(lon_u-lon_d+0.001)
#a = 0
b = lat_u - a*lon_u
#rotate diagonal line
#a = -(lat[0]-lat[2])/(lon[0]-lon[2])
#b = (lat[0] - a*lon[0])

m_lat = (lat[0]+lat[2])/2
m_lon = (lon[0]+lon[2])/2
#plt.scatter(171.4441,-42.7245,c='r') #plotting INZ

df1 = open('list_146_4p8_events_INV.csv')
lines = df1.readlines()
SRF_all = np.zeros([3,len(lines)-1])

df = open('list_146_4p8_events_INV.csv')
for i,x in enumerate(df):
    if(i>0):
        row = x.strip().split(',')
        SRF_all[1,i-1] = row[3]
        SRF_all[2,i-1] = row[2] 
        plt.scatter(float(row[3]),float(row[2]),c='r')
        
        if((float(row[2])-(a*float(row[3])+b))>0 and float(row[2])>lat[0]+0.3):
            print(i)
            print(x)
            SRF_all[0,i-1] = 1
#            SRF.append(row[0])
#            plt.scatter(float(row[3]),float(row[2]),c='r')

#Get the event index with good geometrical distribution:
SRF_index = symmetric_source(SRF_all,m_lon,m_lat,a,b)

df = open('list_146_4p8_events_INV.csv')
SRF_names = []
for i,x in enumerate(df):
    if(i>0) and SRF_index[i-1]==1:
        row = x.strip().split(',')
        if(float(row[2])>lat[0]):
            SRF_names.append(row[0])
            plt.scatter(float(row[3]),float(row[2]),c='g')        

#Save the selected source for inversion:
sname = 'list_'+str(len(SRF_names))+'_4p8_events.txt'
fid=open(sname,'w')
##fid.write("%d\n" %(len(SRF_names)))   
for i in range(0, len(SRF_names)):
    fid.write("%s\n" %(SRF_names[i]))
       
#Save the selected source for inversion with CMT solution:
sname_cmt = 'list_'+str(len(SRF_names))+'_4p8_events.csv'
fid2=open(sname_cmt,'w')

df = open('list_146_4p8_events_INV.csv')
for i,x in enumerate(df):
    if(i==0):
        #First row in the catalog with EventID, Date,Latitude,Longitude,strike1,dip1,rake1, ...
        fid2.write("%s\n" %(x.strip()))
    if(i>0):
        row = x.strip().split(',')
        if(row[0] in SRF_names):
            print(row[0])
            fid2.write("%s\n" %(x.strip()))

for i in range(0,nRec):
    plt.scatter(R[i,0],R[i,1],c='b')   
plt.plot(lon,lat)
#plt.xlim([lonmin,lonmax])
#plt.ylim([latmin,latmax])
plt.xlim([170.5,176])
plt.ylim([-44,-40])
#plt.xlabel('Longitude')
#plt.xlabel('Latitude')
plt.show()
#plt.savefig('VM_BBstations.png',dpi=200)