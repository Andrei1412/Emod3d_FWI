#!/usr/bin/env python
"""
Created on Thu Dec  3 13:40:51 2020

@author: andrei
"""
import numpy as np
import os
import time
from multiprocessing import Pool
import step_length_moduled
from jobs_module import job_submitted
from read_write_module import read_srf_source, read_stat_name, write_par_source_i

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def run_process(process):
    os.system('python {}'.format(process))

# Main script#############
iters = 10

source_file = '../../../StatInfo/SOURCE_SRF.txt'
nShot, S, sNames = read_srf_source(source_file)

station_file = '../../../StatInfo/STATION_dh_4km.txt'
nRec, R, statnames = read_stat_name(station_file)

print('finish creating geometrical correlation')
mkdir_p('../../Vel_es/Vel_es_i')
mkdir_p('../../Vel_ob/Ves_ob_i')
mkdir_p('../../Vel_opt/Vel_ob_i')
mkdir_p('Dump')

ishot_arr = np.linspace(1, nShot, nShot).astype('int')
# ishot_arr=[1]
############################################################
for it in range(1, iters + 1):
    fi = open('iNumber.dat', 'w')
    (np.int64(it)).tofile(fi)
    fi.close()
    print('iNumber=' + str(it))
    #Load the initial model from Model/Models/
    if (it == 1):
        #        os.system('cp ../../../Model/Models/vs3dfile.s ../../../Model/Models/vs3dfile_in.s')
        #        os.system('cp ../../../Model/Models/vp3dfile.p ../../../Model/Models/vp3dfile_in.p')
        #        os.system('cp ../../../Model/Models/rho3dfile.d ../../../Model/Models/rho3dfile_in.d')

        os.system('cp ../../../Model/Models/vs3dfile_2020.s ../../../Model/Models/vs3dfile_in.s')
        os.system('cp ../../../Model/Models/vp3dfile_2020.p ../../../Model/Models/vp3dfile_in.p')
        os.system('cp ../../../Model/Models/rho3dfile_2020.d ../../../Model/Models/rho3dfile_in.d')

    #Load the current model as updated model from the previous iteration in /Dump
    #Make sure the /Dump folder existed to all store data for the current model (Except the synthetic waveform stored in Dump_DD).
    if (it > 1):
        os.system('cp %s' % 'Dump/vs3dfile_iter_' + str(it - 1) + '.s ../../../Model/Models/vs3dfile_in.s')
        os.system('cp %s' % 'Dump/vp3dfile_iter_' + str(it - 1) + '.p ../../../Model/Models/vp3dfile_in.p')
        os.system('cp %s' % 'Dump/rho3dfile_iter_' + str(it - 1) + '.d ../../../Model/Models/rho3dfile_in.d')

    print('Generate synthetic data for pyflex pick:')
    os.system('cp ../../../Model/Models/vs3dfile_in.s Model/vs3dfile_opt.s')
    os.system('cp ../../../Model/Models/vp3dfile_in.p Model/vp3dfile_opt.p')
    os.system('cp ../../../Model/Models/rho3dfile_in.d Model/rho3dfile_opt.d')

    #Remove the kernels and waveform for each individual source from the previous iteration run:
    os.system('rm -r All_shots/*')
    os.system('rm -r ../../Vel_es/*')
    dir_es = '../../Vel_es/Vel_es_i'
    mkdir_p(dir_es)

    #    for ishot in range(1,nShot+1):
    for ishot_id in range(0, len(ishot_arr)):
        ishot = ishot_arr[ishot_id]
        dir_es = '../../Vel_es/Vel_es_' + str(ishot)
        mkdir_p(dir_es)
        dir_opt = '../../Vel_opt/Vel_ob_' + str(ishot)
        mkdir_p(dir_opt)

    job_file = 'fdrun-mpi_Allshot_opt_moduled.sl'
    job_submitted(job_file)
    print('fw emod3d for all sources')

    # input('-->')
    print('Pyflex pick!')
    job_file = 'adj_pyflex_moduled.sl'
    job_submitted(job_file)
    # input('-->back to sum kernels update')
    print('Calculate kernels for all shots')

    #Parallelize the kernel calculations for 2 groups of sources
    #via Master_iterations_part_1_multi_taper.py in PART1/Kernel/Iter/iter1/ e.t.c.
    #In PART1: there're also folders: AdjSims  FwdSims  Kernels for forward/adjoint simulations and kernel calculation.
    for ipart in range(1, 3):

        fi = open('ipart.dat', 'w')
        (np.int64(ipart)).tofile(fi)
        fi.close()
        print('ipart=' + str(ipart))
        #Pass the i_part number (1 or 2) to the corresponding path:
        os.system('mv %s' % 'ipart.dat ../../../PART' + str(ipart) + '/Kernels/Iters/iter1/')

    print('parallel kernel calculation')
    processes = ('Master_iterations_part_1_moduled.py', 'Master_iterations_part_2_moduled.py')
    pool = Pool(processes=2)
    #    processes = ('Master_iterations_part_1_multi_taper.py', 'Master_iterations_part_2_multi_taper.py','Master_iterations_part_3_multi_taper.py', 'Master_iterations_part_4_multi_taper.py')
    #    pool = Pool(processes=4)
    pool.map(run_process, processes)
    #Save the kernels and Hessians as GP_shoti.txt, GS_shoti.txt and HP_shoti.txt, HS_shoti.txt in All_shots:
    print('Finish copy data and kernels')

####Finish parallel calculation of the kernels:
    #Sum kernels and precondition by the approcimated Hessians:
    #Save the summed kernels and Hessians as GP.txt, GS.txt and HP.txt, HS.txt in the running folder (Kernel/Iter/iter1/):
    os.system('python sum_gradients_moduled.py')
    #os.system('python Sum_Gradients_nShot_vsp_save_hessian.py')
    print('finish summing kernels')

    #    #################################
    #L-BFGS gradients calculation (Zhu, 2015).
    #Save gradient in /Dump as Gradp_iter_i.txt and Grads_iter_i.txt
    os.system('python precondition_gradient_L_BFGS_moduled.py')
    print('finish precoditioning kernels')

    ##############################
    #Misfit according to the current model.
    job_file5 = 'observed_misfit.sl'
    job_submitted(job_file5)
    #Save misfit and number of window picked to /Dump as err_iter_i.dat and nwin_iter_i.dat
    f_err0 = open('err_obs.dat', 'r');    Err = np.fromfile(f_err0, dtype=np.float64);    print('Err=' + str(Err))
    os.system('cp err_obs.dat Dump/err_iter_' + str(it) + '.dat')
    f = open('nwin.dat', 'r');    nwin = int(np.fromfile(f, dtype='int64'));    f.close()
    os.system('cp nwin.dat Dump/nwin_iter_' + str(it) + '.dat')
    #Check multi-taper misfit and number of windows picked for the current iteration:
    print('[Err,nwin]=' + str([Err, nwin]))

    #Calculate the optimal step length for updating:
    print('Set the update step length')

    if (it > 2):#Start L-BFGS step-length search if iteration>2.
        #Update model with steplength = 1, save in Model/ as vs3dfile_opt.s, vp3dfile_opt.p and rho3dfile_opt.d:
        job_file = 'model_optimizing_bfgs.sl'
        job_submitted(job_file)
        print('forward simulation with perturbed model')
        #Forward simulation according to the updated model to generate synthetic waveform
        job_file = 'fdrun-mpi_Allshot_opt_moduled.sl'
        job_submitted(job_file)
        print('fw emod3d for all sources')
        #Misfit calculation for the current updated model:
        job_file = 'optimized_misfit.sl'
        job_submitted(job_file)
        print('optimize misfit')
        f_err = open('err_opt.dat', 'r');
        Err2 = np.fromfile(f_err, dtype=np.float64);
        os.system("rm err_opt.dat")
        print('Err2=' + str(Err2))

        #Put constrain on the L-BSGS update: go back to steepest decent method if the update by L-BFGS is too small or misfit does not decrease:
        if (np.abs(Err - Err2) > 0.01 * Err and Err2 < Err):
            # Save the optimal steplength and updated model for the next iteration in /Dump as step_iter_i.dat and
            # vs3dfile_iter_i.s, vp3dfile_iter_i.p and rho3dfile_iter_i.d
            print('Finish update model')
            os.system('mv %s' % 'Model/vs3dfile_opt.s Dump/vs3dfile_iter_' + str(it) + '.s')
            os.system('mv %s' % 'Model/vp3dfile_opt.p Dump/vp3dfile_iter_' + str(it) + '.p')
            os.system('mv %s' % 'Model/rho3dfile_opt.d Dump/rho3dfile_iter_' + str(it) + '.d')

            os.system('mv %s' % 'st.dat Dump/step_iter_' + str(it) + '.dat')
        else:
            #Go back to the stepest decent method if L-BFGS steplength does not sastify the update:
            print('Search for optimal step length')
            opt_step_L2, flag = step_length_moduled.step_length(Err, S, nShot, sNames, it)
            #opt_step_L2, flag = Step_length_moduled.Step_length_new(Err, S, nShot, sNames, it) #for smaller step length.
    else: #Go back to the stepest decent method if iterations<=2.
        print('Search for optimal step length')
        opt_step_L2, flag = step_length_moduled.step_length(Err, S, nShot, sNames, it)
        # opt_step_L2, flag = Step_length_moduled.Step_length_new(Err, S, nShot, sNames, it) #for smaller step length.
    #Save the optimal steplength and updated model for the next iteration in /Dump as step_iter_i.dat and
    #vs3dfile_iter_i.s, vp3dfile_iter_i.p and rho3dfile_iter_i.d
    print('Finish iteration')




