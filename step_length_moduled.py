#!/usr/bin/env python
"""
Created on Thu Dec  3 13:40:51 2020

@author: andrei
"""

import numpy as np
import os
import time
from numpy.linalg import solve
from jobs_module import job_submitted

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def step_length(Err,S,nShot,sNames,it):
    
    L20=Err
    L21=1.1*Err
    L22=0	
    flag=0
    Mu20=0
    Mu21=0
    Mu22=0	
    alt_step=0
    #Initial step length: 2.5% perturbing of the current model:
    st=0.025
#    st=0.05
    count=0
    #Maximum iterations for decreasing/ increasing the step length:
    iter_max=8
          
    while (L21>L20 and count<iter_max and flag==0):
   
        if (count>0):
            st=st/2
        # Update model with decreasing step length, save in Model/ as vs3dfile_opt.s, vp3dfile_opt.p and rho3dfile_opt.d:
        fi=open('st.dat','w'); (np.float64(st)).tofile(fi); fi.close() 
        print('st='+str(st))
        job_file = 'model_optimizing.sl'
        job_submitted(job_file)
        # Forward simulation according to the updated model to generate synthetic waveform
        print('forward simulation with perturbed model')
        job_file = 'fdrun-mpi_Allshot_opt.sl'
        job_submitted(job_file)          
        print('fw emod3d for all sources')
        # Misfit calculation for the current updated model:
        job_file = 'optimized_misfit.sl'
        job_submitted(job_file)          
        print('optimize misfit')
        #Read the current misfit
        f_err=open('err_opt.dat','r'); Err1=np.fromfile(f_err,dtype=np.float64); os.system("rm err_opt.dat")
        L21=Err1
        Mu21=st
        print('L21='+str(L21))
        count=count+1

        if (Mu21<0.01):
            flag=4
            print('flag='+str(flag))
            opt_step_L2=Mu21
            print('st='+str(st))             
            
    if (L21<L20 or Mu21<0.01):
        alt_step=Mu21
    ###
    count=0
    while (L21>L22 and flag==0):
        st=st*1.5
        # Update model with increasing step length, save in Model/ as vs3dfile_opt.s, vp3dfile_opt.p and rho3dfile_opt.d:
        fi=open('st.dat','w'); (np.float64(st)).tofile(fi); fi.close() 
        print('st='+str(st))
        job_file = 'model_optimizing.sl'
        job_submitted(job_file)         
        # Forward simulation according to the updated model to generate synthetic waveform
        print('forward simulation with perturbed model')
        job_file = 'fdrun-mpi_Allshot_opt.sl'
        job_submitted(job_file)          
        print('fw emod3d for all sources')     
        # Misfit calculation for the current updated model:
        job_file = 'optimized_misfit.sl'
        job_submitted(job_file)     
        print('optimize misfit')
        #Read the current misfit
        f_err=open('err_opt.dat','r'); Err1=np.fromfile(f_err,dtype=np.float64); os.system("rm err_opt.dat")
        L22=Err1
        print('L22='+str(L22))
        Mu22=st
        
        if(count>0):
            alt_step=Mu22
        
        count=count+1
        if (count==iter_max) or (Mu22>0.05):

            flag=2
            opt_step_L2=st
        
    if (flag==0):
        A = [[Mu20**2, Mu20, 1], [Mu21**2, Mu21, 1], [Mu22**2, Mu22, 1]]    
        B = [L20,L21,L22]
        print(A)
        print(B)
        X=solve(A,B)
        opt_step_L2 = -X[1]/(2*X[0])
        if (opt_step_L2<0):
            opt_step_L2=Mu20
        if (opt_step_L2>Mu22):
            opt_step_L2=Mu22
        

#Confirm optimal step length
    print('Confirm optimal step length')
    st=opt_step_L2
    # Update model with optimum step length from parabol fitting method, save in Model/ as vs3dfile_opt.s, vp3dfile_opt.p and rho3dfile_opt.d:
    fi=open('st.dat','w'); (np.float64(st)).tofile(fi); fi.close() 
    print('st='+str(st))
    job_file = 'model_optimizing.sl'
    job_submitted(job_file)
    # Forward simulation according to the updated model to generate synthetic waveform
    print('forward simulation with perturbed model')
    job_file = 'fdrun-mpi_Allshot_opt.sl'
    job_submitted(job_file)          
    print('fw emod3d for all sources')
    # Misfit calculation for the current updated model:
    job_file = 'optimized_misfit.sl'
    job_submitted(job_file)
    print('optimize misfit')
    # Read the current misfit
    f_err=open('err_opt.dat','r'); Err2=np.fromfile(f_err,dtype=np.float64); os.system("rm err_opt.dat")
    print('Err2='+str(Err2))

    if (Err2>Err):
        opt_step_L2=alt_step
        print('alt_step='+str(alt_step))
        flag=3
##############################            
    print('flag=')
    print(flag)
    
    print('write updated model')
    st=opt_step_L2
    # Update model with final step length
    fi=open('st.dat','w'); (np.float64(st)).tofile(fi); fi.close() 
    print('st='+str(st))
    job_file = 'model_optimizing.sl'
    job_submitted(job_file)

    #Save the optimal steplength and updated model for the next iteration in /Dump as step_iter_i.dat and
    #vs3dfile_iter_i.s, vp3dfile_iter_i.p and rho3dfile_iter_i.d
    os.system('mv %s' %'Model/vs3dfile_opt.s Dump/vs3dfile_iter_'+str(it)+'.s')
    os.system('mv %s' %'Model/vp3dfile_opt.p Dump/vp3dfile_iter_'+str(it)+'.p')
    os.system('mv %s' %'Model/rho3dfile_opt.d Dump/rho3dfile_iter_'+str(it)+'.d')    
    fi=open('Dump/step_iter_'+str(it)+'.dat','w'); (np.float64(st)).tofile(fi); fi.close()
    return opt_step_L2,flag


def step_length_new(Err, S, nShot, sNames, it):
    #Using smaller step length than the old one.
    L20 = Err
    L21 = 1.1 * Err
    L22 = 0
    flag = 0
    Mu20 = 0
    Mu21 = 0
    Mu22 = 0
    alt_step = 0

    st = 0.025
    count = 0
    iter_max = 8

    while (L21 > L20 and count < iter_max and flag == 0):

        if (count > 0):
            st = st / 2

        fi = open('st.dat', 'w');
        (np.float64(st)).tofile(fi);
        fi.close()
        print('st=' + str(st))
        job_file = 'model_optimizing.sl'
        job_submitted(job_file)

        print('forward simulation with perturbed model')
        job_file = 'fdrun-mpi_Allshot_opt.sl'
        job_submitted(job_file)
        print('fw emod3d for all sources')

        job_file = 'optimized_misfit.sl'
        job_submitted(job_file)
        print('optimize misfit')

        f_err = open('err_opt.dat', 'r');
        Err1 = np.fromfile(f_err, dtype=np.float64);
        os.system("rm err_opt.dat")
        L21 = Err1
        Mu21 = st
        print('L21=' + str(L21))
        count = count + 1

        if (Mu21 < 0.001):
            flag = 4
            print('flag=' + str(flag))
            opt_step_L2 = Mu21
            print('st=' + str(st))

    alt_step = Mu21
    alt_err = Err1
    ###
    count = 0
    while (L21 > L22 and flag == 0):
        st1 = st * 1.5
        fi = open('st.dat', 'w');
        (np.float64(st1)).tofile(fi);
        fi.close()
        print('st=' + str(st))
        job_file = 'model_optimizing.sl'
        job_submitted(job_file)

        print('forward simulation with perturbed model')
        job_file = 'fdrun-mpi_Allshot_opt.sl'
        job_submitted(job_file)
        print('fw emod3d for all sources')

        job_file = 'optimized_misfit.sl'
        job_submitted(job_file)
        print('optimize misfit')

        f_err = open('err_opt.dat', 'r');
        Err1 = np.fromfile(f_err, dtype=np.float64);
        os.system("rm err_opt.dat")
        if (count == 0):
            st = st1
            L22 = Err1
            Mu22 = st

        print('L22=' + str(L22))
        if (Err1 < L22):
            st = st1
            L22 = Err1
            Mu22 = st
        else:
            flag = 1

        count = count + 1
        if (count == iter_max) and (Mu22 > 0.05):
            flag = 2
            opt_step_L2 = st

    if (flag == 0) or (flag == 1):
        A = [[Mu20 ** 2, Mu20, 1], [Mu21 ** 2, Mu21, 1], [Mu22 ** 2, Mu22, 1]]
        B = [L20, L21, L22]
        print(A)
        print(B)
        X = solve(A, B)
        opt_step_L2 = -X[1] / (2 * X[0])
        if (opt_step_L2 < 0):
            opt_step_L2 = Mu20
        if (opt_step_L2 > Mu22):
            opt_step_L2 = Mu22

    # Confirm optimal step length
    print('Confirm optimal step length')
    st = opt_step_L2
    fi = open('st.dat', 'w');
    (np.float64(st)).tofile(fi);
    fi.close()
    print('st=' + str(st))
    job_file = 'model_optimizing.sl'
    job_submitted(job_file)

    print('forward simulation with perturbed model')
    job_file = 'fdrun-mpi_Allshot_opt.sl'
    job_submitted(job_file)
    print('fw emod3d for all sources')

    job_file = 'optimized_misfit.sl'
    job_submitted(job_file)
    print('optimize misfit')

    f_err = open('err_opt.dat', 'r');
    Err2 = np.fromfile(f_err, dtype=np.float64);
    os.system("rm err_opt.dat")
    print('Err2=' + str(Err2))

    if (Err2 > alt_err):
        st = alt_step
        print('alt_step=' + str(alt_step))
        print('write updated model')
        fi = open('st.dat', 'w');
        (np.float64(st)).tofile(fi);
        fi.close()
        print('st=' + str(st))
        job_file = 'model_optimizing.sl'
        job_submitted(job_file)

    print('flag=')
    print(flag)
    os.system('mv %s' % 'Model/vs3dfile_opt.s Dump/vs3dfile_iter_' + str(it) + '.s')
    os.system('mv %s' % 'Model/vp3dfile_opt.p Dump/vp3dfile_iter_' + str(it) + '.p')
    os.system('mv %s' % 'Model/rho3dfile_opt.d Dump/rho3dfile_iter_' + str(it) + '.d')
    fi = open('Dump/step_iter_' + str(it) + '.dat', 'w');
    (np.float64(st)).tofile(fi);
    fi.close()

    return opt_step_L2, flag