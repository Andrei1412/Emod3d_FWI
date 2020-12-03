#!/usr/bin/env python
"""
Created on Thu Dec  3 13:40:51 2020

@author: andrei
"""
import numpy as np
from read_write_module import Vsp_read

def conv_wf_hessian_it(it):
    #emod3d shift for strain tensors
    #flo=0.1
    flo=1.0
    delay_Time=(3/flo)
    ndelay_Time = int(delay_Time/(dt))
    #time delay only for fw wavefield, not bw!
    ef_i=ef[(it+ndelay_Time)*nx*ny*nz:(it+ndelay_Time+1)*nx*ny*nz]
    eb_i=eb[(nt-it-1)*nx*ny*nz:(nt-it)*nx*ny*nz]        
    
    """
         Convolve forward and backward
    """   

    dt_sum=np.zeros(ef_i.shape)
    dt_sum_hessian=np.zeros(ef_i.shape)    
    #slice along y-axis for large chunk of data (out of meamory otherwise!):
    for iy in np.arange(0,ny,1):
        dt_sum[iy*nx*nz:(iy+1)*nx*nz] = np.multiply(ef_i[iy*nx*nz:(iy+1)*nx*nz],eb_i[iy*nx*nz:(iy+1)*nx*nz])  
        dt_sum_hessian[iy*nx*nz:(iy+1)*nx*nz] = np.multiply(eb_i[iy*nx*nz:(iy+1)*nx*nz],eb_i[iy*nx*nz:(iy+1)*nx*nz])         
#    print('it='+str(it))
    return dt_sum, dt_sum_hessian

def conv_wf_hessian(dt):
    #emod3d shift for strain tensors in hessian calculation:
    #flo=0.1
    flo=1.0
    delay_Time=(3/flo)
    ndelay_Time = int(delay_Time/(dt))
    #Sum along time series:
    sum_ij, sum_hessian_ij= conv_wf_hessian_it(0)
    for it in np.arange(1,nt-ndelay_Time,1):
        tmp, tmp_hessian = conv_wf_hessian_it(it)
        sum_ij = sum_ij+dt* tmp
        sum_hessian_ij = sum_hessian_ij+dt* tmp_hessian
        
#    print('done'+eij_name)    
    print([np.min(sum_ij[1]),np.min(sum_hessian_ij[1])])    
    return sum_ij, sum_hessian_ij

def read_strain(eij_name):
    
    matrix_file_fw=fw_file+eij_name
    matrix_file_bw=bw_file+eij_name
    
    fid1=open(matrix_file_fw, 'rb')
    matrix_dummy_fw=np.fromfile(fid1, np.float32)
    matrix_dummy_fw=matrix_dummy_fw[15:]
    
    fid2=open(matrix_file_bw, 'rb')
    matrix_dummy_bw=np.fromfile(fid2, np.float32)  
    matrix_dummy_bw=matrix_dummy_bw[15:]
    
    return matrix_dummy_fw, matrix_dummy_bw

# Read_Dev_strain_ex2.m
nx=88
ny=88
nz=60
nts=1500
nt=300 #strain tensors recorded at every 5 time steps:
#nt=1200
#nt=240

nxyz=nx*ny*nz
#Load strain file parameters:
with open('../../Dev_Strain/fwd01_xyzts.exx', 'rb') as fid:
    data_array = np.fromfile(fid, np.float32)
    data_array=data_array[0:14]
    dx=data_array[8]
    dy=data_array[9]
    dz=data_array[10]
    dt=data_array[11] #dt=5 * time sampling
    del data_array

print('dt='+str(dt))    
#Load input structure model:
snap_file1='../../../../Model/Models/vs3dfile_in.s'
snap_file2='../../../../Model/Models/vp3dfile_in.p'
snap_file3='../../../../Model/Models/rho3dfile_in.d'

Vs=Vsp_read(nx,ny,nz,snap_file1)
Vp=Vsp_read(nx,ny,nz,snap_file2)
rho=Vsp_read(nx,ny,nz,snap_file3)
#Converst to Lame's parameters:
Lamda=np.multiply(rho,(np.multiply(Vp,Vp)-2*np.multiply(Vs,Vs)));
Mu=np.multiply(rho,np.multiply(Vs,Vs))
Kappa=Lamda+2/3*Mu

#6 strain components:
eijp_names = {'exx','eyy','ezz','exy','exz','eyz'};
fw_file='../../Dev_Strain/fwd01_xyzts.';
bw_file='../../Dev_Strain/adj01_xyzts.';

for eij_name in eijp_names:
    #Read strain components:
    print(eij_name)
    if eij_name=='exx':
        exx3, exxb = read_strain(eij_name)
        
    if eij_name=='eyy':
        eyy3, eyyb = read_strain(eij_name)

    if eij_name=='ezz':
        ezz3, ezzb = read_strain(eij_name)
        
    if eij_name=='exy':
        ef,eb = read_strain(eij_name)
        tmp1, tmp1_hessian = conv_wf_hessian(dt)
        sum_xy=4*tmp1
        sum_hessian_xy=4*tmp1_hessian        
        
    if eij_name=='exz':
        ef,eb = read_strain(eij_name)
        tmp1, tmp1_hessian = conv_wf_hessian(dt)
        sum_xz=4*tmp1
        sum_hessian_xz=4*tmp1_hessian                    
        
    if eij_name=='eyz':
        ef,eb = read_strain(eij_name)
        tmp1, tmp1_hessian = conv_wf_hessian(dt)
        sum_yz=4*tmp1
        sum_hessian_yz=4*tmp1_hessian      

#sum of 3 normal strains:
ef = exx3+eyy3+ezz3
eb = exxb+eyyb+ezzb
tmp1, tmp1_hessian = conv_wf_hessian(dt)
sum1 = tmp1
sum1_hessian = tmp1_hessian
#Kappa kernel:
Kappa_arr=Kappa.reshape(-1)
GK =  np.multiply(sum1,Kappa_arr)
HK =  np.multiply(sum1_hessian,Kappa_arr)

#Mu kernel:
sum2 = 0.5*(sum_xy+sum_xz+sum_yz)  
sum2_hessian = 0.5*(sum_hessian_xy+sum_hessian_xz+sum_hessian_yz)      
Mu_arr=Mu.reshape(-1)    

GM = 2*np.multiply(sum2,Mu_arr)
HM = 2*np.multiply(sum2_hessian,Mu_arr)

#Convert to kernels and hessians for Vs and Vp:
GS=2*(GM-4/3*np.multiply(np.divide(Mu_arr,Kappa_arr),GK))
GP=2*np.multiply(np.divide((Kappa_arr+4/3*Mu_arr),Kappa_arr),GK)	  

Hessian_S=2*(HM-4/3*np.multiply(np.divide(Mu_arr,Kappa_arr),HK))
Hessian_P=2*np.multiply(np.divide((Kappa_arr+4/3*Mu_arr),Kappa_arr),HK)	      
#Save kernels and hessians for Vs and Vp:
np.savetxt('KS.txt', GS)
np.savetxt('KP.txt', GP)

np.savetxt('HS.txt', Hessian_S)
np.savetxt('HP.txt', Hessian_P)
