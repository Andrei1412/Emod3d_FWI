#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 09:56:20 2019

@author: user
"""
import os
from subprocess import Popen, PIPE
import numpy as np
#import matplotlib.pyplot as plt
def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def read_srf(srf_file):
    """
    Function for reading srf files
    """
    with open(srf_file, 'r') as f:
        lines = f.readlines()
        line5=lines[5].split()
        Si=np.zeros((1,3))
        Si[0,0]=line5[0]
        Si[0,1]=line5[1]
        Si[0,2]=line5[2]
    
    return Si

def read_list_srf(list_file):

    with open(list_file, 'r') as f:
        lines = f.readlines()
    nSource = len(lines) 
    sNames = [] 
    for i in range(0,nSource,1):
        line_i=lines[i].split('\n')
        sNames.append(line_i[0])
        
    return nSource, sNames


def ll2xy_conv(MODEL_LON,MODEL_LAT,MODEL_ROT,ELON,ELAT,EDEP,NX,NY,HH):
       
    XAZIM = MODEL_ROT+90.0
    
    cmd = 'echo '+ str(ELON) + " " + str(ELAT) +'| ./ll2xy mlon='+str(MODEL_LON)+' mlat='+str(MODEL_LAT)+' xazim='+str(XAZIM)
    stdout = Popen(cmd, shell=True, stdout=PIPE).stdout
    output = stdout.read()    
    EXY = output.split()
   # input('-->')

    XSRC = int(0.5*NX + float(EXY[0])/HH)

    YSRC = int(0.5*NY + float(EXY[1])/HH)

    ZSRC = int(EDEP/HH + 0.5) + 1 
    
    print([XSRC, YSRC, ZSRC])   
    
    return XSRC, YSRC, ZSRC

def write_srf_source(S,sNames,S_LL):
    filename1='SOURCE_SRF.txt' 
    Ns=len(S)
    print(filename1)
    fid = open(filename1,'w')
    fid.write("%4d\n" %(Ns))    
    count=0
    for i in range(0,Ns,1):
        fid.write("%4d%4d%4d %20s %20f%20f%10f\n" %(S[count,0],S[count,1],S[count,2],sNames[count],S_LL[count,0],S_LL[count,1],S_LL[count,2]))        
        count=count+1    
    

#########################################
nSource, sNames = read_list_srf('list_13_4p8_events.txt')
        
S = np.zeros((nSource,3))
S_LL = np.zeros((nSource,3))

#Given from VM defined for emod3d run:
NX = 88
NY = 88
NZ = 60
HH = 4.0

MODEL_LON = 173.099917317
MODEL_LAT = -42.0999452931
MODEL_ROT = 0.0

#Find all srf-source format for 13 events in the list.
#Otherwise generate the srf-file from Geonet catalog (CMT solution given in csv-file of the same name, i.e "list_13_4p8_events.txt.csv")
for i in range(0,nSource,1):
    dir_srf = 'SRF_'+str(nSource)+'s'
    mkdir_p(dir_srf)
    #Andrei_Sources_20191209: folder contains all srf-files for events > Mw 4 in the domain created by Robin.
    os.system('cp Andrei_Sources_20191209/'+sNames[i]+'/Srf/'+sNames[i]+'.srf '+dir_srf)
    print('cp Andrei_Sources_20191209/'+sNames[i]+'/Srf/'+sNames[i]+'.srf '+dir_srf)
    srf_file = dir_srf+'/'+sNames[i]+'.srf'    
    
    print(srf_file)
    #input('-->')
    try:
        #Read the lon/lat/depth from srf file for the according event and convert to Cartesian coordinates:
        Si = read_srf(srf_file)
        ELON = Si[0,0]
        ELAT = Si[0,1]
        EDEP = Si[0,2]
        
        S_LL[i,:] = [ELON, ELAT, EDEP]
    
        S[i,:] = ll2xy_conv(MODEL_LON,MODEL_LAT,MODEL_ROT,ELON,ELAT,EDEP,NX,NY,HH)
    except:
        print('no srf found in 422 events')
    
write_srf_source(S,sNames,S_LL)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
