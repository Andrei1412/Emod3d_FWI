#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
#import os
import numpy as np
import os
import time

from jobs_module import job_submitted
from read_write_module import read_srf_source, write_par_opt_source_i

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


#Run for the first itertion only:
iters=1
source_file='../../../StatInfo/SOURCE_SRF.txt'
#source_file='../../../StatInfo/SOURCE_checker.txt'
nShot, S, sNames = read_srf_source(source_file)

mkdir_p('FWD/All_par')
mkdir_p('FWD/All_srf')

############################################################
E=np.zeros((iters,1))
for it in range(1,iters+1):

    print('iNumber='+str(it))
    #Model files for forward simulation only:
    if (it==1):
        os.system('cp ../../../Model/Models/vs3dfile.s Model/vs3dfile_opt.s')
        os.system('cp ../../../Model/Models/vp3dfile.p Model/vp3dfile_opt.p')
        os.system('cp ../../../Model/Models/rho3dfile.d Model/rho3dfile_opt.d')

    if (it>1):
        os.system('cp %s' %'Dump/vs3dfile_iter_'+str(it-1)+'.s Model/vs3dfile_opt.s')
        os.system('cp %s' %'Dump/vp3dfile_iter_'+str(it-1)+'.p Model/vp3dfile_opt.p')
        os.system('cp %s' %'Dump/rho3dfile_iter_'+str(it-1)+'.d Model/rho3dfile_opt.d')
        
    time.sleep(5)
    #Generate par-file and srf-file for
    ishot_arr=np.linspace(1,nShot,nShot).astype('int')

    for ishot_id in range(0,len(ishot_arr)):
        ishot=ishot_arr[ishot_id]

        print('isource='+str(ishot))
        #Edit e3d.par file with source location:
        write_par_opt_source_i(S[ishot-1,:])
        os.system('cp FWD/e3d_mysource_i_opt.par FWD/All_par/e3d_mysource_'+str(ishot)+'_opt.par')
        #Store srf file for each source:
        os.system('cp SRF_13s/'+str(sNames[ishot-1])+'.srf FWD/All_srf/srf_source_'+str(ishot)+'.srf')

        dir_opt='../../Vel_opt/Vel_ob_'+str(ishot)
        mkdir_p(dir_opt)

    input('check all_srf')

    #run this slurm for multi-event forward modeling.
    # To check for log data, set set SOURCES = ( 1 ), run fw simulation and check at Rlog/
    job_file = 'fdrun-mpi_Allshot_opt_moduled.sl'
    job_submitted(job_file)  

#Check dir_opt for waveform output:
     
    

    
