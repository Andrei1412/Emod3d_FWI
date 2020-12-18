#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:00:22 2019
Plot the observed/ synthetic waveforms and pyflex/ pyadjoint results
@author: andrei
"""
import obspy
import pyflex
import pprint
import logging
import pandas as pd
from scipy.stats import norm
from obspy.io.sac import SACTrace
#from pysac.sactrace import SACTrace
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import fftpack
from scipy import signal
from qcore import timeseries
from scipy import integrate
from scipy.signal import butter, lfilter
import pyadjoint
#from pyadjoint import multitaper_misfit as mt
#import multitaper_misfit_new

from subprocess import Popen, PIPE

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def winpad(lt,t_off,t_on,pad):
    """
    Sin taper window
    """    
    #pad=5
    L=t_off-t_on+2*pad
    
#    window=signal.gaussian(30, 4)
#    window1=signal.resample(window,L,axis=0, window=None)
#    window=np.linspace(1, 1, L)))
    window=np.ones((L))
    
    #x=np.arange(0,pad,1)
    x=np.linspace(0, np.pi/2, pad)
    sinx=np.sin(x)
    window[0:pad] = sinx
    window[L-pad:L] = sinx[::-1]    
    print('lt='+str(lt))    
    ar1=np.zeros((t_on-pad))
    ar2=np.zeros((lt-t_off-pad))
    window_pad0 = np.concatenate((ar1,window))
    window_pad = np.concatenate((window_pad0,ar2))    
    
    return window_pad  
############################################################test_sac_flexwin_148events.py

def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def readGP_2(loc, fname):
    """
    Convinience function for reading files in the Graves and Pitarka format
    """
    with open("/".join([loc, fname]), 'r') as f:
        lines = f.readlines()
    
    data = []

    for line in lines[2:]:
        data.append([float(val) for val in line.split()])

    data=np.concatenate(data) 
    
    line1=lines[1].split()
    num_pts=int(line1[0])
    dt=float(line1[1])
    shift=float(line1[4])

    return data, num_pts, dt, shift

def read_stat_name_utm(station_file):

    with open(station_file, 'r') as f:
        lines = f.readlines()
    line0=lines[0].split()
    nRec=int(line0[0])
    R=np.zeros((nRec,2))
    statnames = [] 
    for i in range(1,nRec+1):
        line_i=lines[i].split()
        R[i-1,0]=float(line_i[3])
        R[i-1,1]=float(line_i[4])
        statnames.append(line_i[2])
    return nRec, R, statnames

def read_source_new_utm(source_file):
    
    with open(source_file, 'r') as f:
        lines = f.readlines()    
    line0=lines[0].split()
    nShot=int(line0[0])
    S=np.zeros((nShot,3))
    sNames=[]
    for i in range(1,nShot+1):
        line_i=lines[i].split()
        S[i-1,0]=float(line_i[4])
        S[i-1,1]=float(line_i[5])
        S[i-1,2]=float(line_i[6])
        sNames.append(line_i[3])
    
    return nShot, S, sNames

def rms(stat_data):

    num_pts=len(stat_data)
    D = (np.sum(np.square(stat_data))/num_pts)**0.5

    return stat_data/D

def time_shift(data,t_shift,dt):
    n_pts = len(data)
    nshift_T = int(t_shift/(dt))
    data_shift = np.zeros(data.shape)
    data_shift[nshift_T:n_pts] = data[0:n_pts-nshift_T]
    return data_shift

def time_shift_emod3d(data,delay_Time,dt):
    n_pts = len(data)
    ndelay_Time = int(delay_Time/(dt))
    data_shift = np.zeros(data.shape)
    data_shift[0:n_pts-ndelay_Time] = data[ndelay_Time:n_pts]
    return data_shift
    
def write_win_qual(R_Time_record_end,dt,windows,filename):
    fid = open(filename,'w')
    fid.write("%s\n" % ('# NUM_WIN =')) 
    fid.write("%s\n" % ('# i win_start win_end Tshift CC dlnA') )
 
    for i in range(0,len(windows)):
        t_off=windows[i].right*dt 
        if(R_Time_record_end<windows[i].right*dt):
            t_off=R_Time_record_end
        fid.write("%4d%20f%20f%20f%20f%20f\n" %(i+1,windows[i].left*dt,t_off,windows[i].cc_shift*dt,windows[i].max_cc_value,windows[i].dlnA))
    return    

def write_adj_source_ts(v1, mainfolder_source, source, dt):
    #write adjoint source/ ascii file of name A.x e.g.
    # filename1=mainfolder_source+v1
    vs1 = v1.split('.')
    timeseries.seis2txt(source, dt, mainfolder_source, vs1[0], vs1[1])
    return
    
###############################################################
station_file = 'STATION_utm_dh_4km.txt'
nRec, R, statnames = read_stat_name_utm(station_file)
source_file = 'SOURCE_SRF.txt'
nShot, S, sNames = read_source_new_utm(source_file)

num_pts=1500; dt=0.16;
#num_pts=3000; dt=0.08;
t = np.arange(num_pts)*dt

############/nesi/nobackup/nesi00213/RunFolder/tdn27/rgraves/Adjoint/Syn_VMs/Kernels/#########################
fs = 1/dt
#lowcut = 0.05
highcut = 0.1
lowcut = 0.025
#highcut = 0.075

R_Time_record = np.zeros([2,nShot,nRec])
R_Time_record[0,:,:]=0;R_Time_record[1,:,:]=240;

#flo=0.1
flo=1.0
delay_Time=(3/flo)
t_shift=0
t_min = 0;
t_max = max(t);
#
GV=['.090','.000','.ver']
BH=['BHX','BHY','BHZ']

GV=['.090','.000','.ver']
GV_ascii=['.x','.y','.z']

mainfolder='/home/andrei/workspace/GMPlots/Sim/Vel_es/*'
mainfolder_o='/home/andrei/workspace/GMPlots/Sim/Vel_ob/*'

#Uncomment this to see all window picked parameters
logger = logging.getLogger("pyflex")
logger.setLevel(logging.DEBUG)
# DEFAULT pyflex value for Marlborough data set
config = pyflex.Config(
    min_period=1 / highcut, max_period=1 / lowcut,
    # min_period=10.0, max_period=40.0,
    stalta_waterlevel=0.08, tshift_acceptance_level=10.0,
    #    stalta_waterlevel=0.08, tshift_acceptance_level=5.0,
    dlna_acceptance_level=5.0, dlna_reference=0.0, cc_acceptance_level=0.7,
    #    dlna_acceptance_level=5.0, dlna_reference=0.0,cc_acceptance_level=0.9,

    s2n_limit=3.0,
    #    min_surface_wave_velocity=3.0, max_time_before_first_arrival=-10.0,
    c_0=0.7, c_1=2.0, c_2=0.0, c_3a=3.0, c_3b=2.0, c_4a=2.5, c_4b=12.0,
    check_global_data_quality=True, snr_integrate_base=5.0, snr_max_base=3.5,
    noise_start_index=0, noise_end_index=None,
    signal_start_index=None, signal_end_index=-1, window_weight_fct=None,
    window_signal_to_noise_type=u'amplitude', resolution_strategy=u'interval_scheduling')

# DEFAULT pyadjoint value for Marlborough data set
config_adj = pyadjoint.Config(
    measure_type='dt',  # travel time adj source
    #        measure_type='am',#amplitude abnormaly adj source
    #        measure_type='wf',#full waveform adj source
    min_period=1 / highcut, max_period=1 / lowcut,
    # min_period=10.0, max_period=40.0,
    lnpt=15,
    transfunc_waterlevel=1.0E-10,
    water_threshold=0.02,
    ipower_costaper=10,
    min_cycle_in_window=0.5,
    taper_type='hann',
    taper_percentage=0.1,
    mt_nw=4.0,
    num_taper=5,
    dt_fac=2.0,
    phase_step=1.5,
    err_fac=2.5,
    dt_max_scale=3.5,
    #        measure_type='dt',
    dt_sigma_min=1.0,
    dlna_sigma_min=0.5,
    use_cc_error=False,
    use_mt_error=False)

filename = 'WIN.SAC.win.qual'
#for ishot in range(1,nShot+1):
#for ishot in range(1,1+1):
td_max_count=[] 
index_all=np.zeros((nRec,3,nShot)) 

#ishot_arr=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
ishot_arr=[1]
os.system('rm ALL_WINs_pyflex_temp/*')
os.system('rm -r /home/andrei/workspace/GMPlots/Sim/Vel_MT_adj/*')

for ishot_id in range(0,len(ishot_arr)):
    ishot=ishot_arr[ishot_id]
    print('ishot='+str(ishot))
    
    mkdir_p('/home/andrei/workspace/GMPlots/Sim/Vel_MT_adj/Vel_ob_'+str(ishot))
    mainfolder_source='/home/andrei/workspace/GMPlots/Sim/Vel_MT_adj/Vel_ob_'+str(ishot)+'/'
    mainfolder='/home/andrei/workspace/GMPlots/Sim/Vel_es/'
           
    os.system('rm /home/andrei/workspace/GMPlots/Sim/Vel_ob/*')
    os.system('rm /home/andrei/workspace/GMPlots/Sim/Vel_es/*')    
    
    os.system('scp maui:/scale_wlg_nobackup/filesets/nobackup/nesi00213/RunFolder/tdn27/rgraves/NZVMs/Marlborough_Events_4KM/NewVM_20200207/INV_MarlVM_DH_4km_template/Kernels/Vel_ob_240s_HH_13s/Vel_ob_'+str(ishot)+'/*.* /home/andrei/workspace/GMPlots/Sim/Vel_ob/')     
    os.system('scp maui:/scale_wlg_nobackup/filesets/nobackup/nesi00213/RunFolder/tdn27/rgraves/NZVMs/Marlborough_Events_4KM/NewVM_20200207/INV_MarlVM_DH_4km_template/Kernels/Vel_es/Vel_es_'+str(ishot)+'/*.* /home/andrei/workspace/GMPlots/Sim/Vel_es/')  

    #input('-->')
    for i,statname in enumerate(statnames):
        plt.figure(figsize=(15,1.25))  
        for k in range(0,3):

            # input time serie s0/ adjoint time serie v0:
            s0 = statname + GV[k]
            v0 = statname + GV_ascii[k]
            # Read synthetic data for the current model and correct by emod3d shift:
            data_sim = timeseries.read_ascii('/home/andrei/workspace/GMPlots/Sim/Vel_es/' + s0)
            data_sim = time_shift_emod3d(data_sim, delay_Time, dt)

            # Assign adjoint source as zero time serie:
            adj_source = np.zeros(np.shape(data_sim))

            try:
                data_obs = timeseries.read_ascii('/home/andrei/workspace/GMPlots/Sim/Vel_ob/' + s0)
            except:
                # write zero-adjoint with an exception in observed data reading:
                write_adj_source_ts(v0, mainfolder_source, adj_source, dt)
                continue  # continue here if no record for this station/ component

            # Tapering:
            data_obs = np.multiply(signal.tukey(int(num_pts), 0.05), data_obs)
            data_sim = np.multiply(signal.tukey(int(num_pts), 0.05), data_sim)
            # Filtering
            data_obs = butter_bandpass_filter(data_obs, lowcut, highcut, fs, order=4)
            data_sim = butter_bandpass_filter(data_sim, lowcut, highcut, fs, order=4)
            # Integrating:
            data_obs = np.cumsum(data_obs) * dt
            data_sim = np.cumsum(data_sim) * dt

            # Formating to sac:
            obs_sac = SACTrace(kstnm=statname, kcmpnm=BH[k], stla=R[i, 1], stlo=R[i, 0], evla=S[ishot - 1, 1],
                               evlo=S[ishot - 1, 0], evdp=S[ishot - 1, 2], nzyear=2000, nzjday=1, nzhour=0, nzmin=0,
                               nzsec=0, nzmsec=0,
                               t0=t_min, t1=t_max, delta=dt, b=t_shift, data=data_obs)
            obs_sac.write('OBS.SAC', byteorder='little')

            syn_sac = SACTrace(kstnm=statname, kcmpnm=BH[k], stla=R[i, 1], stlo=R[i, 0], evla=S[ishot - 1, 1],
                               evlo=S[ishot - 1, 0], evdp=S[ishot - 1, 2], nzyear=2000, nzjday=1, nzhour=0, nzmin=0,
                               nzsec=0, nzmsec=0,
                               t0=t_min, t1=t_max, delta=dt, b=t_shift, data=data_sim)
            syn_sac.write('SYN.SAC', byteorder='little') 
            
            plt.figure(figsize=(10,2.5))    
            plt.plot(t,data_obs,c='k',linestyle='solid')
            plt.plot(t,data_sim,c='r',linestyle='solid')
#            plt.ylim([-0.01, 0.01])
            plt.show()
            print(s0)
#                    
            obs_data = obspy.read("OBS.SAC")
            synth_data = obspy.read("SYN.SAC")
            
            obs_data.detrend("linear")
            obs_data.taper(max_percentage=0.05, type="hann")
#            obs_data.filter("bandpass", freqmin=1.0 / 40.0, freqmax=10.0 / 2.0,
#                            corners=4, zerophase=True)
            
            synth_data.detrend("linear")
            synth_data.taper(max_percentage=0.05, type="hann")
#            synth_data.filter("bandpass", freqmin=1.0 / 40.0, freqmax=10.0 / 2.0,
#                              corners=4, zerophase=True)                    
#            
            #Select windows
            try:
#                    windows = pyflex.select_windows(obs_data, synth_data, config, plot=True)  
                windows = pyflex.select_windows(obs_data, synth_data, config, plot=True)     
#                 windows = flexwin_new.select_windows(obs_data1, synth_data1, config, plot=False)                     

            except:
                windows = []
                            
#            adj_src = calculate_adjoint_source(obs_data[0].data, synth_data[0].data, dt, config_adj, windows,adjoint_src = True, figure = False)
            
            if (len(windows) > 0):
                #                if(1>0):
                # Remove window less than 30s:
                if ((windows[0].right - windows[0].left) * dt > 30):
                    write_win_qual(R_Time_record[1,ishot-1,i],dt,windows,filename)
                    print(len(windows))
                    index_all[i,k,ishot-1] = 1
                
                    try:
                        config_adj.measure_type = 'dt'                           
                        adj_src = pyadjoint.calculate_adjoint_source("multitaper_misfit",obs_data,synth_data,config_adj, windows,adjoint_src=True,plot=True)
#                        if(adj_src.misfit<50):
                        if(1>0):
                            adj_source = adj_src.adjoint_source   
                        
                        print('total misfit:'+str(adj_src.misfit))
                        
#                        config_adj.measure_type = 'am'                        
#                        adj_src = pyadjoint.calculate_adjoint_source("multitaper_misfit",obs_data,synth_data,config_adj, windows,adjoint_src=True,plot=True)
#                        adj_source = adj_src.adjoint_source     +          adj_source      
#                        print('total misfit:'+str(adj_src.misfit))                        
                    
                    except:
                        print('raise Exception')     
        ##                        
#                plt.plot(adj_source);plt.show();  
                write_adj_source_ts(v0, mainfolder_source, adj_source, dt)
#                        
#                e_s_c_name = str(ishot)+'.'+s0+'.win'
#                os.system('cp WIN.SAC.win.qual ALL_WINs_pyflex_temp/'+e_s_c_name)      


