#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
import os
from subprocess import Popen, PIPE
import numpy as np
import math
import matplotlib.pyplot as plt
#import pdb; 
def xy2ll_conv(MODEL_LON,MODEL_LAT,MODEL_ROT,XSRC,YSRC,NX,NY,HH):
       
    XAZIM = MODEL_ROT+90.0
    
    EXY0 = (XSRC-0.5*NX)*HH
    EXY1 = (YSRC-0.5*NY)*HH
    
    cmd = 'echo '+ str(EXY0) + " " + str(EXY1) +'| ./xy2ll mlon='+str(MODEL_LON)+' mlat='+str(MODEL_LAT)+' xazim='+str(XAZIM)
    stdout = Popen(cmd, shell=True, stdout=PIPE).stdout
    output = stdout.read()    
    ELL = output.split()
   # input('-->')
    ELON = ELL[0]
    ELAT = ELL[1]    
   
    print([ELON, ELAT])   
    
    return ELON, ELAT

def ll2xy_conv(MODEL_LON,MODEL_LAT,MODEL_ROT,ELON,ELAT,NX,NY,HH):
       
    XAZIM = MODEL_ROT+90.0
       
    cmd = 'echo '+ str(ELON) + " " + str(ELAT) +'| ./ll2xy mlon='+str(MODEL_LON)+' mlat='+str(MODEL_LAT)+' xazim='+str(XAZIM)
    stdout = Popen(cmd, shell=True, stdout=PIPE).stdout
    output = stdout.read()    
    EXY = output.split()
   # input('-->')

    XSRC = int(0.5*NX + float(EXY[0])/HH)

    YSRC = int(0.5*NY + float(EXY[1])/HH)
   
    print([XSRC, YSRC])   
    
    return XSRC, YSRC


##########################################################
NX = 88
NY = 88
NZ = 60
HH = 4.0

MODEL_LON = 173.099917317
MODEL_LAT = -42.0999452931
MODEL_ROT = 0.0

nx=88
ny=88
nz=60
dx=4; dy=4;

INZ_lon = 171.4441
INZ_lat = -42.7245
#INZ_lon = 171.437546
#INZ_lat = -42.698776

[INZ_x, INZ_y] = ll2xy_conv(MODEL_LON,MODEL_LAT,MODEL_ROT,INZ_lon,INZ_lat,NX,NY,HH)
[INZ_lon_check, INZ_lat_check] = xy2ll_conv(MODEL_LON,MODEL_LAT,MODEL_ROT,INZ_x,INZ_y,NX,NY,HH)

print([INZ_x, INZ_y] )
print([INZ_lon_check, INZ_lat_check] )

