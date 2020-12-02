#!/usr/bin/env python
# Sum kernels and precondition by the approcimated Hessians:
#Save the summed kernels and Hessians as GP.txt, GS.txt and HP.txt, HS.txt
import numpy as np
import scipy.ndimage as ndimage
from read_write_module import read_srf_source

def threshole_hessian(Hessian, THRESHOLD_HESS):
    [ny, nz, nx] = Hessian.shape
    Inv_Hessian = np.zeros([ny, nz, nx])
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                if (Hessian[j, k, i] > THRESHOLD_HESS):
                    Inv_Hessian[j, k, i] = 1 / Hessian[j, k, i]
                else:
                    Inv_Hessian[j, k, i] = 1 / THRESHOLD_HESS

    return Inv_Hessian

def tape_xyz_matrix(nx,ny,nz,R_xyz):
    #Tape boundaries
    Tp=np.ones((ny,nz,nx))
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            for k in range(1,nz+1):

                if (i<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*i/R_xyz)

                if ((nx-i-1)<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*(nx-i-1)/R_xyz)

                if (j<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*j/R_xyz)

                if ((ny-j-1)<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*(ny-j-1)/R_xyz)
                #Tape near bottom boundary
                if ((nz-k-1)<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*(nz-k-1)/R_xyz)

    return Tp

#Model dimensions:
nx = 88
ny = 88
nz = 60

#Crearte taper matrix for R_nf grids from left. right, east, west and bottom and save to /Dump as Tp_xyz.txt:
R_nf = 5
Tp=tape_xyz_matrix(nx,ny,nz,R_nf);Tp_arr=Tp.reshape(-1);
tape_file='Dump/Tp_xyz.txt';np.savetxt(tape_file, Tp_arr)

source_file = '../../../StatInfo/SOURCE_SRF.txt'
nShot, S, _ = read_srf_source(source_file)
##########################################################
#Gradient matrices
Gra_S = np.zeros((ny, nz, nx))
Gra_P = np.zeros((ny, nz, nx))
#Hessian matrices
Hessian_S = np.zeros((ny, nz, nx))
Hessian_P = np.zeros((ny, nz, nx))

#Source array:
ishot_arr = np.linspace(1, nShot, nShot).astype('int')
#Gather kernels and hessians for all events to calculate the gradients:
for ishot_id in range(0, len(ishot_arr)):
    #load the kernels for each source:
    ii = ishot_arr[ishot_id]

    print('Sum Gradient Shot ' + str(ii))
    GS_file = 'All_shots/GS_shot' + str(ii) + '.txt'
    GP_file = 'All_shots/GP_shot' + str(ii) + '.txt'

    GS_arr = np.loadtxt(GS_file)
    GS = np.reshape(GS_arr, [ny, nz, nx])
    GP_arr = np.loadtxt(GP_file)
    GP = np.reshape(GP_arr, [ny, nz, nx])

    # load the hessians for each source:
    print('Sum Hessian Shot ' + str(ii))
    HS_file = 'All_shots/HS_shot' + str(ii) + '.txt'
    HP_file = 'All_shots/HP_shot' + str(ii) + '.txt'

    HS_arr = np.loadtxt(HS_file)
    HS = np.reshape(HS_arr, [ny, nz, nx])
    HP_arr = np.loadtxt(HP_file)
    HP = np.reshape(HP_arr, [ny, nz, nx])

    #No tapering near source/ receiver as Romanovics, 1996:
    Tp = np.ones((ny, nz, nx))
    #Gradients as sum of the kernels
    Gra_S = Gra_S + np.multiply(GS, Tp)
    Gra_P = Gra_P + np.multiply(GP, Tp)
    #Hessian preconditioners as sum of the absolute hessians for each source:
    Hessian_S = Hessian_S + np.abs(np.multiply(HS, Tp))
    Hessian_P = Hessian_P + np.abs(np.multiply(HP, Tp))

##Tape boundary: 5-grid from the boundaries:
tape_file = 'Dump/Tp_xyz.txt'
Tp_arr = np.loadtxt(tape_file)
Tp = np.reshape(Tp_arr, [ny, nz, nx])
#
Gra_S = np.multiply(Gra_S, Tp)
Gra_P = np.multiply(Gra_P, Tp)

Hessian_S = np.multiply(Hessian_S, Tp)
Hessian_P = np.multiply(Hessian_P, Tp)

#Nullify the gradients at the surface to account for the non-linear in the waveform at the surface:
Gra_S[:, 0, :] = 0
Gra_P[:, 0, :] = 0

Hessian_S[:, 0, :] = 0
Hessian_P[:, 0, :] = 0
#Normalize the hessian preconditioners:
Hessian_S = Hessian_S / np.max(Hessian_S)
Hessian_P = Hessian_P / np.max(Hessian_P)

# Smooth the gradients/ hessians before invert (Specfem3d manual)
#Spatial Gaussian function with radius of 8x4x8 km in x,z,y directions:
Gra_S = -ndimage.gaussian_filter(Gra_S, [2, 1, 2])
Gra_P = -ndimage.gaussian_filter(Gra_P, [2, 1, 2])

Hessian_S = ndimage.gaussian_filter(Hessian_S, [2, 1, 2])
Hessian_P = ndimage.gaussian_filter(Hessian_P, [2, 1, 2])
#Threshhold for matrix inversion:
THRESHOLD_HESS = 5.e-4
Inv_Hessian_S = threshole_hessian(Hessian_S, THRESHOLD_HESS)
Inv_Hessian_P = threshole_hessian(Hessian_P, THRESHOLD_HESS)

#Precondition the gradients with hessian preconditioners:
Gra_S = np.multiply(Gra_S, Inv_Hessian_S)
Gra_P = np.multiply(Gra_P, Inv_Hessian_P)

#Save the summed gradient (steepest gradient for model update/ CG or L-BFGS implementation)
Gra_S_arr = Gra_S.reshape(-1)
np.savetxt('Gra_S.txt', Gra_S_arr)

Gra_P_arr = Gra_P.reshape(-1)
np.savetxt('Gra_P.txt', Gra_P_arr)
print('finish Summing Gradients')

print('max Gs=' + str(np.max(Gra_S_arr)))
print('max Gp=' + str(np.max(Gra_P_arr)))

