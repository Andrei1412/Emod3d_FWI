#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 11:24:47 2019

@author: user
"""
import numpy as np
from numpy import linalg as LA
import os
from read_write_module import Vsp_read, write_model_Vsp

#
nx=88
ny=88
nz=60

#Read iteration number
fi1=open('iNumber.dat','r')
it=int(np.fromfile(fi1,dtype='int64'))
fi1.close()
print('iNumber='+str(it))

#Read steplength for l-bfgs
fi1=open('st.dat','r')
st=np.fromfile(fi1,dtype=np.float64)
fi1.close()
print('st='+str(st))

#Load current models for i-th iteration:
snap_file1='../../../Model/Models/vs3dfile_in.s'
snap_file2='../../../Model/Models/vp3dfile_in.p'
snap_file3='../../../Model/Models/rho3dfile_in.d'
Vs=Vsp_read(nx,ny,nz,snap_file1)
Vp=Vsp_read(nx,ny,nz,snap_file2)
rho=Vsp_read(nx,ny,nz,snap_file3)

#Load current gradients for i-th iteration:
print('Dump/Grads_iter_'+str(it)+'.txt')    
Gra_S_arr = np.loadtxt('Dump/Grads_iter_'+str(it)+'.txt')   
Gra_P_arr = np.loadtxt('Dump/Gradp_iter_'+str(it)+'.txt')  

Vsp=np.zeros(2*nx*ny*nz)
grad_Vsp=np.zeros(2*nx*ny*nz)

Vsp[0:nx*ny*nz] = Vs.reshape(-1)
Vsp[nx*ny*nz:2*nx*ny*nz] = Vp.reshape(-1)

grad_Vsp[0:nx*ny*nz] = Gra_S_arr
grad_Vsp[nx*ny*nz:2*nx*ny*nz] = Gra_P_arr

#Update models Vs/Vp and constrain the update if needed:
Vsp0=Vsp-st*grad_Vsp

Vp_max=9.0
Vp_min=1.0
Vs_max=4.5
Vs_min=0.5

for i in range(1,nx*ny*nz+1):

    if Vsp0[nx*ny*nz+i-1]>Vp_max:
        Vsp0[nx*ny*nz+i-1]=Vp_max
    if Vsp0[nx*ny*nz+i-1]<Vp_min:
        Vsp0[nx*ny*nz+i-1]=Vp_min

    if Vsp0[i-1]>Vs_max:
        Vsp0[i-1]=Vs_max
    if Vsp0[i-1]<Vs_min:
        Vsp0[i-1]=Vs_min

Vs0=np.reshape(Vsp0[0:nx*ny*nz],[ny,nz,nx])
Vp0=np.reshape(Vsp0[nx*ny*nz:2*nx*ny*nz],[ny,nz,nx])

print('Mu21='+str(st))	
print('max st_Vp0*Gra_Vp ='+str(np.max(np.abs(st*grad_Vsp))))
#Save the updated models for forward modeling at Model/ as vs3dfile_opt.s, vp3dfile_opt.p and  rho3dfile_opt.d
write_model_Vsp(rho,Vs0,Vp0)
#time.sleep(5)
os.system('mv vs3dfile1.s Model/vs3dfile_opt.s')
os.system('mv vp3dfile1.p Model/vp3dfile_opt.p')
os.system('mv rho3dfile1.d Model/rho3dfile_opt.d')
