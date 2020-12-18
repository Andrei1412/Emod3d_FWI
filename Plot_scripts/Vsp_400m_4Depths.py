# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:54:45 2019

@author: user
"""

#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
#import os
import numpy as np
import scipy.ndimage as ndimage
from matplotlib import cm
import matplotlib.pyplot as plt

#import pdb; 
def Vsp_read(nx,ny,nz,snap_file):
    fid=open(snap_file,'r')
    sek=np.fromfile(fid, dtype='<f4')
    Vp=np.reshape(sek,[ny,nz,nx])
    return Vp   

# #############################
nx=880
ny=880
nz=300
#nz=120

dx=0.4
#X = np.linspace(0,nx-1) * dx
#Y = np.linspace(0,ny-1) * dx
x_positions  = np.linspace(0,nx-1,5)
x_labels = np.linspace(0,nx-1,5) * dx
y_positions = np.linspace(0,ny-1,5)
y_labels  = np.linspace(0,ny-1,5) * dx

#snap_file1='vs3dfile.s'
#snap_file2='vp3dfile.p'
#snap_file3='rho3dfile.d'

#snap_file1='vs3dfile_int_4_2km.s'
#snap_file2='vp3dfile_int_4_2km.p'
#snap_file3='rho3dfile_int_4_2km.d'

snap_file1='vs3dfile_inv_4km_400m.s'
snap_file2='vp3dfile_inv_4km_400m.p'
snap_file3='rho3dfile_inv_4km_400m.d'


Vs1=Vsp_read(nx,ny,nz,snap_file1)
Vp1=Vsp_read(nx,ny,nz,snap_file2)
rho1=Vsp_read(nx,ny,nz,snap_file3)

#snap_file1='vs3dfile_iter_5.s'
#snap_file2='vp3dfile_iter_5.p'
#snap_file3='rho3dfile_iter_5.d'
#
#Vs=Vsp_read(nx,ny,nz,snap_file1)
#Vp=Vsp_read(nx,ny,nz,snap_file2)
#rho=Vsp_read(nx,ny,nz,snap_file3)

#Vs1_Vs0 = np.true_divide(Vs,Vs1)
#fid10=open('Vs1_Vs0.s','wb')
#sek2 = np.array(Vs, dtype=np.float32)
#sek2.astype('float32').tofile(fid10)
z0_arr = [21, 41, 61, 91]
#z0_arr = [1, 2, 3, 4]




fig = plt.figure(figsize=(20,16))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
m = cm.ScalarMappable(cmap=cm.jet)
#m.set_array([-0.2,0.2])    
#m.set_array([-0.1,0.1]) 
for i in range(0, len(z0_arr)):
    z0 = z0_arr[i]
    
#    LL = np.true_divide(Vs[:,z0,:],Vs1[:,z0,:])
#    LL = np.true_divide(Vs,Vs1)
    #LL3D = np.true_divide(Vs,Vs1)
    LL_log = Vs1[:,z0,:]
    
#    LL_log = np.log(LL)
#    LL_log=ndimage.gaussian_filter(LL_log, 5)   
#    LL_log=ndimage.gaussian_filter(LL_log, 10)       

    ax = fig.add_subplot(2, 2, i+1)
#    ax.text(0.5, 0.5, str('z='+str((z0-1)*4)+'km'),fontsize=14, ha='center')

    #plt.style.use('classic')
    plt.imshow(LL_log,cmap = cm.jet_r)
    plt.xticks(x_positions, x_labels)
    plt.yticks(y_positions, y_labels)    
    plt.xlabel('x (km)',fontsize=20,)
    plt.ylabel('y (km)',fontsize=20,)    
#    plt.imshow(X,Y,LL_log)    
#    im = ax.imshow(LL_log)
#    plt.colorbar(im,ax)    
    plt.title('z='+str((z0-1)*dx)+'km',fontsize=20,)
#    plt.clim(-0.1,0.1)
#    plt.clim(-0.2,0.2)    
    
    plt.colorbar().set_label('Vs [km/s]', labelpad=-10, y=1.2, rotation=0,fontsize=20,)
#    plt.colorbar(m).set_label('log(Vs1/Vs)', labelpad=-10, y=1.2, rotation=0)
#plt.show()
plt.savefig('4depths_Vs.png',dpi=100)

