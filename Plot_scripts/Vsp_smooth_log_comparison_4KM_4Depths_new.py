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

#snap_file1='3366146/vs3dfile.s'
#snap_file2='3366146/vp3dfile.p'
#snap_file3='3366146/rho3dfile.d'

#snap_file1='vs3dfile_smooth_4km.s'
#snap_file2='vp3dfile_smooth_4km.p'
#snap_file3='rho3dfile_smooth_4km.d'
#
snap_file1='NZVM_2020/vs3dfile_2020.s'
snap_file2='NZVM_2020/vp3dfile_2020.p'
snap_file3='NZVM_2020/rho3dfile_2020.d'

#filepath1 = 'NZVM_2020/NZVM2020_Dump_29_09_20_15iters_13s_BB_data_G212_bfgs_dt_adj_ntape5_dtts5_cc07_good/'
#snap_file1=filepath1+'vs3dfile_iter_9.s'
#snap_file2=filepath1+'vp3dfile_iter_9.p'
#snap_file3=filepath1+'rho3dfile_iter_9.d'

Vs1=Vsp_read(nx,ny,nz,snap_file1)
Vp1=Vsp_read(nx,ny,nz,snap_file2)
rho1=Vsp_read(nx,ny,nz,snap_file3)

snap_file1='vs3dfile_iter_5.s'
snap_file2='vp3dfile_iter_5.p'
snap_file3='rho3dfile_iter_5.d'
#snap_file1='../Dump_AmbientNoise/AB_Model_10BB_15_11_20/vs3dfile_iter_5.s'
#snap_file2='../Dump_AmbientNoise/AB_Model_10BB_15_11_20/vp3dfile_iter_5.p'
#snap_file3='../Dump_AmbientNoise/AB_Model_10BB_15_11_20/rho3dfile_iter_5.d'

#snap_file1='NZVM_2020/vs3dfile_2020.s'
#snap_file2='NZVM_2020/vp3dfile_2020.p'
#snap_file3='NZVM_2020/rho3dfile_2020.d'

#filepath1 = 'NZVM_2020/NZVM2020_Dump_29_09_20_15iters_13s_BB_data_G212_bfgs_dt_adj_ntape5_dtts5_cc07_good/'
#filepath1 = '../Dump_AmbientNoise/Checker_10BB_test/'

#snap_file1=filepath1+'vs3dfile_iter_3.s'
#snap_file2=filepath1+'vp3dfile_iter_3.p'
#snap_file3=filepath1+'rho3dfile_iter_3.d'

#snap_file1=filepath1+'vs3dfile_true_square_4x_G5x4.s'
#snap_file2=filepath1+'vp3dfile_true_square_4x_G5x4.p'
#snap_file3=filepath1+'rho3dfile_true_square_4x_G5x4.d'

#snap_file1='../RELOCATION_EVENTS_INV/vs3dfile_iter_2.s'
#snap_file2='../RELOCATION_EVENTS_INV/vp3dfile_iter_2.p'
#snap_file3='../RELOCATION_EVENTS_INV/rho3dfile_iter_2.d'

#file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/Dump_01_10_20_3iters_27s_BB_data_G212_bfgs_dt_cc09_Tmax_5s_NZVM2020_2nd_run_good/'
#file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/Dump_04_10_20_5iters_27s_BB_data_G212_bfgs_dt_cc08_Tmax_5s_NZVM2020_M10_2nd_run/'
#file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/Dump_04_10_20_AB_SYN/'
#file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/Dump_02_10_20_2iters_27s_BB_data_G212_bfgs_dt_cc09_Tmax_5s_NZVM2020_M9_2nd_run_good/'
#file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/NZVM2010_Dump_29_09_20_15iters_13s_BB_data_G212_bfgs_dt_adj_ntape5_dtts5_cc07_good_NZVM_2020/'
##file_path='/home/andrei/workspace/GMPlots/Marlborough_Events_4KM/INVERSION_GOOD_NZVM2020/'
#snap_file1=file_path+'vs3dfile_iter_1.s'
#snap_file2=file_path+'vp3dfile_iter_1.p'
#snap_file3=file_path+'rho3dfile_iter_1.d'

Vs=Vsp_read(nx,ny,nz,snap_file1)
Vp=Vsp_read(nx,ny,nz,snap_file2)
rho=Vsp_read(nx,ny,nz,snap_file3)

#Vs1_Vs0 = np.true_divide(Vs,Vs1)
#fid10=open('Vs1_Vs0.s','wb')
#sek2 = np.array(Vs, dtype=np.float32)
#sek2.astype('float32').tofile(fid10)
z0_arr = [2, 3, 4, 5]
#z0_arr = [1, 2, 3, 4]




fig = plt.figure(figsize=(20,16))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
#m = cm.ScalarMappable(cmap=cm.jet)
m = cm.ScalarMappable(cmap='bwr')
m.set_array([-0.2,0.2])    
#m.set_array([-0.1,0.1]) 
for i in range(0, len(z0_arr)):
    z0 = z0_arr[i]
    
    LL = np.true_divide(Vs[:,z0,:],Vs1[:,z0,:])
    #LL = np.true_divide(Vs,Vs1)
    #LL3D = np.true_divide(Vs,Vs1)
    
    LL_log = np.log(LL)
#    LL_log=ndimage.gaussian_filter(LL_log, 3)   
#    LL_log=ndimage.gaussian_filter(LL_log, 10)       

    ax = fig.add_subplot(2, 2, i+1)
#    ax.text(0.5, 0.5, str('z='+str((z0-1)*4)+'km'),fontsize=14, ha='center')

    #plt.style.use('classic')
#    plt.imshow(LL_log,cmap = cm.jet_r)
    plt.imshow(LL_log,cmap = 'bwr_r')    
    plt.xticks(x_positions, x_labels)
    plt.yticks(y_positions, y_labels)    
    plt.xlabel('x (km)',fontsize=20,)
    plt.ylabel('y (km)',fontsize=20,)    
#    plt.imshow(X,Y,LL_log)    
#    im = ax.imshow(LL_log)
#    plt.colorbar(im,ax)    
    plt.title('z='+str((z0-1)*dx)+'km',fontsize=20,)
    plt.clim(-0.2,0.2)
#    plt.clim(-0.05,0.05)    
    
    plt.colorbar().set_label('log(Vs1/Vs)', labelpad=-10, y=1.2, rotation=0,fontsize=20,)
#    plt.colorbar(m).set_label('log(Vs1/Vs)', labelpad=-10, y=1.2, rotation=0)
#plt.show()
plt.savefig('4depths_Vs1_Vs.png',dpi=300)

