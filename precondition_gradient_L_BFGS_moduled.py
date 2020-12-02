#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:52:35 2019

@author: user
"""
import numpy as np
from numpy import linalg as LA
from read_write_module import Vsp_read

def get_gradient(istore):
    #istore = current iteration
    #load current gradient matrices and vectorize them for computation efficiency
    grads_istore = np.loadtxt('Dump/Grads_iter_' + str(istore) + '.txt')
    gradp_istore = np.loadtxt('Dump/Gradp_iter_' + str(istore) + '.txt')
    gradient = np.zeros(2 * nx * ny * nz)
    gradient[0:nx * ny * nz] = grads_istore
    gradient[nx * ny * nz:2 * nx * ny * nz] = gradp_istore

    return gradient

def get_model(istore):
    #load current velocity models and vectorize them for computation efficiency
    #No density update required
    snap_file1 = 'Dump/vs3dfile_iter_' + str(istore) + '.s'
    snap_file2 = 'Dump/vp3dfile_iter_' + str(istore) + '.p'
    Vs = Vsp_read(nx, ny, nz, snap_file1)
    Vp = Vsp_read(nx, ny, nz, snap_file2)
    vs_istore = Vs.reshape(-1)
    vp_istore = Vp.reshape(-1)

    model = np.zeros(2 * nx * ny * nz)
    model[0:nx * ny * nz] = vs_istore
    model[nx * ny * nz:2 * nx * ny * nz] = vp_istore

    return model


def get_model_init():
    #load initial velocity models and vectorize them for computation efficiency
    #    snap_file1='../../../Model/Models/vs3dfile.s'
    #    snap_file2='../../../Model/Models/vp3dfile.p'
    snap_file1 = '../../../Model/Models/vs3dfile_2020.s'
    snap_file2 = '../../../Model/Models/vp3dfile_2020.p'

    Vs = Vsp_read(nx, ny, nz, snap_file1)
    Vp = Vsp_read(nx, ny, nz, snap_file2)
    vs_istore = Vs.reshape(-1)
    vp_istore = Vp.reshape(-1)

    model = np.zeros(2 * nx * ny * nz)
    model[0:nx * ny * nz] = vs_istore
    model[nx * ny * nz:2 * nx * ny * nz] = vp_istore

    return model


###########################################################################
nx = 88
ny = 88
nz = 60

#Load the current iteration number from Master loop:
fi1 = open('iNumber.dat', 'r')
it = int(np.fromfile(fi1, dtype='int64'))
fi1.close()
print('iNumber=' + str(it))

#Load the steepest gradient and vectorize them for computation efficiency
Gra_S_arr = np.loadtxt('Gra_S.txt')
Gra_Vs = np.reshape(Gra_S_arr, [ny, nz, nx])

Gra_P_arr = np.loadtxt('Gra_P.txt')
Gra_Vp = np.reshape(Gra_P_arr, [ny, nz, nx])

gra = np.zeros(2 * nx * ny * nz)
gra[0:nx * ny * nz] = Gra_S_arr
gra[nx * ny * nz:2 * nx * ny * nz] = Gra_P_arr

# Preconditioning/ no-Hessian kernels
#Zhu, 2015; Nocedal, 2006 or specfem3d manual for more detail of the formulations:
# gra = gra/LA.norm(gra)

if (it > 1):
    #Load model from one iteration back:
    model_old = get_model(it - 1)
else:
    model_old = get_model_init()
#  ! initialize arrays
#  a(:)=0.0
#  p(:)=0.0
#  gradient1(:)=0.0
#  gradient0(:)=0.0
#  model1(:)=0.0
#  model0(:)=0.0
#  gradient_diff(:)=0.0
#  model_diff(:)=0.0
#  q_vector(:)=0.0
#  r_vector(:)=0.0
a = np.zeros(it)
p = np.zeros(it)

# Current gradient at it-iteration
#q_vector as current gradient
q_vector = gra
#Choose number of iterations for l-bfgs, min(m)=3
# m_store = 3
m_store = 5
# *******starting backward store *****************'

# gra1=gra
if (it <= 2):
    #r_vector as steppest gradients for two first iterations:
    r_vector = gra

if (it > 2):
    #load C, and their derivatives, inverted value of the summed product
    iter_store = it - m_store
    if (iter_store < 1):
        iter_store = 1

    for istore in range(it - 1, iter_store, -1):
        gradient1 = get_gradient(int(istore))
        gradient0 = get_gradient(int(istore - 1))

        model1 = get_model(int(istore))
        model0 = get_model(int(istore - 1))

        gradient_diff = gradient1 - gradient0
        model_diff = model1 - model0
        #inverted value of the summed product between gradient difference and model difference:
        p_sum = np.sum(np.multiply(gradient_diff, model_diff))
        p[istore] = 1.0 / p_sum
        #the summed product  of q_vector and model difference; multiply by p-value
        a_sum = np.sum(np.multiply(q_vector, model_diff))
        a[istore] = p[istore] * a_sum
        print('a,p,istore:' + str([a[istore], p[istore], istore]))
        #Update q_vector with a-value from reversed count
        q_vector = q_vector - (gradient_diff * a[istore])

    #load the last iteration's gradients, models
    istore = it - 1
    gradient1 = get_gradient(int(istore))
    gradient0 = get_gradient(int(istore - 1))

    model1 = get_model(int(istore))
    model0 = get_model(int(istore - 1))

    gradient_diff = gradient1 - gradient0
    model_diff = model1 - model0

    # ! this implements Algorithm 7.4 and equation (7.20) on page 178 of the book of
    # ! Jorge Nocedal and Stephen Wright, "Numerical Optimization", Springer, second edition (2006)
    p_k_up = np.sum(np.multiply(gradient_diff, model_diff))
    p_k_down = np.sum(np.multiply(gradient_diff, gradient_diff))
    print('p_k_up,p_k_down:' + str([p_k_up, p_k_down]))

    #Update p_k and r_vector for it>2
    p_k = p_k_up / p_k_down
    r_vector = p_k * q_vector

    # ********starting forward store ***********
    #Load q_vector
    for istore in range(iter_store + 1, it):
        gradient1 = get_gradient(int(istore))
        gradient0 = get_gradient(int(istore - 1))

        model1 = get_model(int(istore))
        model0 = get_model(int(istore - 1))

        gradient_diff = gradient1 - gradient0
        model_diff = model1 - model0

        b_sum = np.sum(np.multiply(gradient_diff, r_vector))
        b = p[istore] * b_sum
        # Update r_vector with b-value from forward count
        r_vector = r_vector + model_diff * (a[istore] - b)

        print('b,istore:' + str([b, istore]))
    #Negative direction applied:
    r_vector = -1.0 * r_vector

#################################################################
#Reshape the l-bfgs gradients and save to Dump/ as Gradp_iter_i.txt, Grads_iter_i.txt
Gra_Vs_arr = r_vector[0:nx * ny * nz]
Gra_Vp_arr = r_vector[nx * ny * nz:2 * nx * ny * nz]

print('write gradient iter' + str(it))
np.savetxt('Dump/Grads_iter_' + str(it) + '.txt', Gra_Vs_arr)
np.savetxt('Dump/Gradp_iter_' + str(it) + '.txt', Gra_Vp_arr)

#L-BFGS step length:
epsilon = 0.01
if (it <= 2):
    #1% perturbing of the model if it<=2 (ignore this as optimized step length was calculated out of this function using stepest decent direction)
    st = epsilon * LA.norm(model_old) / LA.norm(r_vector)

if (it > 2):
    #Default step length = 1 for l-bfgs gradients:
    st = 1
    #    st = np.max([1,epsilon*LA.norm(model_old)/LA.norm(r_vector)])
    print([1, epsilon * LA.norm(model_old) / LA.norm(r_vector)])
#Save the L-BFGS step length:
fi = open('st.dat', 'w');
(np.float64(st)).tofile(fi);
fi.close()
print('st=' + str(st))
