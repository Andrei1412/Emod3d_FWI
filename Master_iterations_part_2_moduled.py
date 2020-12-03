#!/usr/bin/env python
"""
Created on Thu Dec  3 13:40:51 2020

@author: andrei
"""
import os
import numpy as np
import time
from jobs_module import job_submitted
from read_write_module import read_srf_source, read_stat_name, write_par_source_i

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

# Save the main running path as Kernels/Iters/iter1/:
pwd = os.getcwd()
print(os.getcwd())
# Go to main folder:
os.chdir('../../../')
path_main = os.getcwd()
# Go to PART1/Kernels/Iters/iter1/ as running folder:
#path = '../../../PART1/Kernels/Iters/iter1/'
path = path_main+'/PART2/Kernels/Iters/iter1/'
os.chdir(path)
print(os.getcwd())

#Read source file in cartesian coordinates
source_file = path_main+'/StatInfo/SOURCE_SRF.txt'
nShot, S, sNames = read_srf_source(source_file)
#Read station file in cartesian coordinates
station_file = path_main+'/StatInfo/STATION_dh_4km.txt'
nRec, R, statnames = read_stat_name(station_file)

#Prepare folder for simmulated waveform from the current
mkdir_p('../../Vel_es/Vel_es_i')
############################################################
#Load the index matrix of channels with picked windows:
R_all_arr=np.loadtxt(path_main+'/Kernels/index_all_ncc_pyflex.txt')
R_all=R_all_arr.reshape([nRec,3,nShot])
#Load the ipart (1 or 2) for the current running folder.
fi1=open('ipart.dat','r')
ipart=np.int64(np.fromfile(fi1,dtype='int64'))
fi1.close()
ipart=int(ipart)

#Move the script from main running folder Kernels/Iters/iter1/ to the current running folder PARTi/Kernels/Iters/iter1/
os.system('cp '+path_main+'/Kernels/Iters/iter1/FWT_emod3d_shot_i_part1.sl FWT_emod3d_shot_i_part1.sl')
os.system('cp '+path_main+'/Kernels/Iters/iter1/e3d_mysource_xyz_default.par.part e3d_mysource_xyz_default.par')
os.system('cp '+path_main+'/Kernels/Iters/iter1/set_run+merge_params_new_h.csh.part ADJ/set_run+merge_params_new_h.csh')
os.system('cp '+path_main+'/Kernels/Iters/iter1/FWT_emod3d_shot_i_part2.sl FWT_emod3d_shot_i_part2.sl')
os.system('cp '+path_main+'/Kernels/Iters/iter1/ascii2adj_ker_iter.csh ../../../AdjSims/V3.0.7-a2a_xyz/ascii2adj_ker_iter.csh')
os.system('cp '+path_main+'/Kernels/Iters/iter1/kernel_shot_i_hessian.sl kernel_shot_i_hessian.sl')
os.system('cp '+path_main+'/Kernels/Iters/iter1/calc_kernels_moduled.py calc_kernels_vsp_hessian.py')

#Kernels calculation of a group of events:
ishot_arr=[7,8,9,10,11,12,13]

for ishot_id in range(0,len(ishot_arr)):
    ishot=ishot_arr[ishot_id]

    if ((np.sum(R_all[:,:,ishot-1])>0)):
        # Check if the number of window picked for this event is larger than 0
        # (the actual pick should be larger than 10 for 30 channels from 10 broadband stations).
        fi=open('iShot.dat','w')
        (np.int64(ishot)).tofile(fi)
        fi.close() 
        print('isource='+str(ishot))
        #Clean data for current iteration according to event i in the main path:
        os.system('rm -r '+path_main+'/Kernels/Vel_es/Vel_es_'+str(ishot)+'/*')
        os.system('rm -r '+path_main+'/Kernels/Iters/iter1/All_shots/GS_shot'+str(ishot)+'.txt')
        os.system('rm -r '+path_main+'/Kernels/Iters/iter1/All_shots/GP_shot'+str(ishot)+'.txt')
        os.system('rm ../../Dev_Strain/*')
        os.system('rm -r ../../../AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii/*')
        #Write e3d.par file from template FWD/e3d_mysource_xyz_default.par for Emod3d: FwdSims/V3.0.7-xyz (forward simulation with waveform and forward strain tensor output)
        write_par_source_i(S[ishot-1,:])
        #Assign the srf-file for the current source from  Kernels/Iters/iter1/SRF_13s
        os.system('cp '+path_main+'/Kernels/Iters/iter1/SRF_13s/'+str(sNames[ishot-1])+'.srf srf_file.srf')
        #Run Emod3d: FwdSims/V3.0.7-xyz (forward simulation with waveform and forward strain tensor output) for a single event:
        job_file11 = 'FWT_emod3d_shot_i_part1.sl'
        job_submitted(job_file11)
        #Copy adjoint source calculated in advance from  Kernels/Vel_MT_adj to AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii/ for adjoint simulation
        os.system('cp '+path_main+'/Kernels/Vel_MT_adj/Vel_ob_'+str(ishot)+'/*.* ../../../AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii/')
        print("adjoint source calculation finished")
    
        #Run Emod3d: AdjSims/V3.0.7-a2a_xyz (adjoint simulation with backward strain tensor output) for a single event:
        #e3d.par file for this simulation was generated by a csh-script: ADJ/set_run+merge_params_new_h.csh
        job_file12 = 'FWT_emod3d_shot_i_part2.sl'
        job_submitted(job_file12)

        #Calculate kernels and estimated Hessians for the current events using python  script calc_kernels_vsp_hessian.py:
        job_file2 = 'kernel_shot_i_hessian.sl'
        job_submitted(job_file2)
    #    print("kernel calculation finished")
        time.sleep(5)
        #Save the kernels to Kernels/Iters/iter1/All_shots and waveform to Kernels/Vel_es in the main path:
        os.system('mv ../../Vel_es/Vel_es_i/*.* '+path_main+'/Kernels/Vel_es/Vel_es_'+str(ishot))
        os.system('mv KS.txt '+path_main+'/Kernels/Iters/iter1/All_shots/GS_shot'+str(ishot)+'.txt')
        os.system('mv KP.txt '+path_main+'/Kernels/Iters/iter1/All_shots/GP_shot'+str(ishot)+'.txt')
        os.system('mv HS.txt '+path_main+'/Kernels/Iters/iter1/All_shots/HS_shot'+str(ishot)+'.txt')
        os.system('mv HP.txt '+path_main+'/Kernels/Iters/iter1/All_shots/HP_shot'+str(ishot)+'.txt')
        
    else: #If not enough picked channels for the current event then no simullation/ kernel calcultion was performed.
        print('No kernel calculated for source '+str(ishot))             
#Go back to the main running path as Kernels/Iters/iter1/ before finish kernel calculation for PART1:
os.chdir(pwd)
print(os.getcwd())    
