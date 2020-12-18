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
nx=88
ny=88
nz=60

dx=4
#X = np.linspace(0,nx-1) * dx
#Y = np.linspace(0,ny-1) * dx
x_positions  = np.linspace(0,nx-1,5)
x_labels = np.linspace(0,nx-1,5) * dx
y_positions = np.linspace(0,ny-1,5)
y_labels  = np.linspace(0,ny-1,5) * dx

#snap_file1='vs3dfile_smooth_4km.s'
#snap_file2='vp3dfile_smooth_4km.p'
#snap_file3='rho3dfile_smooth_4km.d'

#snap_file1='NZVM_2020/vs3dfile_2020.s'
#snap_file2='NZVM_2020/vp3dfile_2020.p'
#snap_file3='NZVM_2020/rho3dfile_2020.d'

#snap_file1='NZVM_2020/vs3dfile_iter_10.s'
#snap_file2='NZVM_2020/vp3dfile_iter_10.p'
#snap_file3='NZVM_2020/rho3dfile_iter_10.d'

#file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/Dump_01_10_20_3iters_27s_BB_data_G212_bfgs_dt_cc09_Tmax_5s_NZVM2020_2nd_run_good/'
#file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/Dump_04_10_20_5iters_27s_BB_data_G212_bfgs_dt_cc08_Tmax_5s_NZVM2020_M10_2nd_run/'
file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/Dump_04_10_20_AB_SYN/'#
snap_file1=file_path+'vs3dfile_iter_1.s'
snap_file2=file_path+'vp3dfile_iter_1.p'
snap_file3=file_path+'rho3dfile_iter_1.d'
#
Vs=Vsp_read(nx,ny,nz,snap_file1)
Vp=Vsp_read(nx,ny,nz,snap_file2)
rho=Vsp_read(nx,ny,nz,snap_file3)

#Vs1_Vs0 = np.true_divide(Vs,Vs1)
#fid10=open('Vs1_Vs0.s','wb')
#sek2 = np.array(Vs, dtype=np.float32)
#sek2.astype('float32').tofile(fid10)
z0_arr = [1, 2, 4, 6]
#fig = plt.figure(figsize=(10,10))

Vs_mean_arr = [3.16055,3.42018,3.6864,3.98812]


fig = plt.figure(figsize=(20,16))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
m = cm.ScalarMappable(cmap=cm.jet)
#m.set_array([-0.2,0.2])    
m.set_array([-0.1,0.1]) 
for i in range(0, len(z0_arr)):
    z0 = z0_arr[i]
    
#    LL = np.true_divide(Vs[:,z0,:],Vs1[:,z0,:])
    #LL = np.true_divide(Vs,Vs1)
    #LL3D = np.true_divide(Vs,Vs1)
    
#    LL_log = np.log(LL)
    LL_log = Vs[:,z0,:]

    Vs_mean=np.mean(LL_log)
#    Vs_mean = Vs_mean_arr[i]
    print('Vs_mean='+str(Vs_mean))
    
#    LL_log=ndimage.gaussian_filter(LL_log, 5)   

    ax = fig.add_subplot(2, 2, i+1)
#    ax.text(0.5, 0.5, str('z='+str((z0-1)*4)+'km'),fontsize=14, ha='center')

    #plt.style.use('classic')
    plt.imshow(LL_log,cmap = cm.jet_r)
    plt.xticks(x_positions, x_labels)
    plt.yticks(y_positions, y_labels)    
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')    
#    plt.imshow(X,Y,LL_log)    
#    im = ax.imshow(LL_log)
#    plt.colorbar(im,ax)    
    plt.title('z='+str((z0-1)*4)+'km')
#    plt.clim(-0.1,0.1)
#    plt.clim(Vs_mean*0.8,Vs_mean*1.2)    
    
    plt.colorbar().set_label('Vs [km/s]', labelpad=-10, y=1.2, rotation=0, fontsize=14, ha='center')
#    plt.colorbar(m).set_label('log(Vs1/Vs)', labelpad=-10, y=1.2, rotation=0)
#plt.show()
plt.savefig('4depths_Vs.png',dpi=100)

