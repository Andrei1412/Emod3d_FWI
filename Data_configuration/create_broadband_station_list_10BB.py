#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:33:28 2019

@author: user
"""
import numpy as np
import os

def read_stat_name(station_file):

    with open(station_file, 'r') as f:
        lines = f.readlines()
    line0=lines[0].split()
    nRec=int(line0[0])
    R=np.zeros((nRec,3))
    statnames = [] 
    for i in range(1,nRec+1):
        line_i=lines[i].split()
        R[i-1,0]=int(line_i[0])
        R[i-1,1]=int(line_i[1])
        R[i-1,2]=int(line_i[2])
        statnames.append(line_i[3])
    return nRec, R, statnames

def read_list_station_ll(list_file):

    with open(list_file, 'r') as f:
        lines = f.readlines()
    nStat = len(lines)
    stat_names = []
    R_list=np.zeros((nStat,2))

    for i in range(0,nStat,1):
        line_i=lines[i].split()

        R_list[i,0]=float(line_i[0])
        R_list[i,1]=float(line_i[1])
        stat_names.append(line_i[2])

    return nStat, R_list, stat_names

def write_BB_station(statnames_BB,R_cartesian):
    filename1='STATION_dh_4km.txt'
    Nr=len(R_cartesian)
    print(filename1)
    fid = open(filename1,'w')
    fid.write("%4d\n" %(Nr))
    count=0
    for i in range(0,Nr,1):
        fid.write("%4d%4d%4d%20s\n" %(R_cartesian[count,0],R_cartesian[count,1],1, statnames_BB[count]))
        count=count+1

    return

def write_utm_BB_station(statnames_BB,R_LL):
    filename1='STATION_utm_dh_4km.txt'
    Nr=len(R_LL)
    print(filename1)
    fid = open(filename1,'w')
    fid.write("%4d\n" %(Nr))
    count=0
    for i in range(0,Nr,1):
        fid.write("%4d%4d %20s %20f%20f\n" %(R_LL[count,0],R_LL[count,1],statnames_BB[count],R_LL[count,2],R_LL[count,3]))
        count=count+1

    return

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
        
###########################
#List of 10 Broadband stations located in the domain:
statnames_BB = ['QRZ','NNZ','WEL','DSZ','THZ','KHZ','INZ','LTZ','GVZ','OXZ']        
#statnames_BB = ['QRZ','NNZ','WEL','DSZ','THZ','KHZ','LTZ','GVZ','OXZ']          
        
nRec_BB = len(statnames_BB)
R_BB = np.zeros((nRec_BB,4))

#Write cartesian coord. for BB stations from list of all strong motion stations:
station_file = 'fd_rt01-h4.000.statcords'    
_, R, statnames = read_stat_name(station_file)

for i,stat_BB in enumerate(statnames_BB):
    if(stat_BB in statnames):
        x = statnames.index(stat_BB)
        R_BB[i,0] = R[x,0]
        R_BB[i,1] = R[x,1]
    #Since no INZ station listed in the strong motion station list, an adhoc script created for this one only
    #lon/lat for INZ station can be found at: https://www.geonet.org.nz/data/network/sensor/INZ
    if(stat_BB=='INZ'): #from Check_INZ.py
        R_BB[i,0] = 10
        R_BB[i,1] = 61      
        
R_cartesian = R_BB[:,0:2]
#Write BB stations in Cartesian only
write_BB_station(statnames_BB,R_cartesian)

#Write utm coordinates (lon/lat) for BB stations from list of all strong motion stations:
station_file_ll = 'fd_rt01-h4.000.ll'    
_, R_list, stat_names = read_list_station_ll(station_file_ll)

for i,stat_BB in enumerate(statnames_BB):
    if(stat_BB in stat_names):
        x = stat_names.index(stat_BB)
        R_BB[i,2] = R_list[x,0]
        R_BB[i,3] = R_list[x,1]
        
    if(stat_BB=='INZ'): #from Check_INZ.py
        R_BB[i,2] = 171.4441
        R_BB[i,3] = -42.7245              
#Write BB stations in Cartesian and lon/lat
write_utm_BB_station(statnames_BB,R_BB)        


           
