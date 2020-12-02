#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 14:19:45 2018

@author: andrei
"""
import obspy
import pyflex
import pyadjoint
import flexwin_new
from obspy.io.sac import SACTrace

import numpy as np
import os
from qcore import timeseries
from scipy import signal

from read_write_module import read_source_new_utm, read_stat_name_utm
from time_serie_module import winpad, butter_bandpass_filter, rms, time_shift, time_shift_emod3d, write_win_qual, \
    write_adj_source_ts



###################################################
# READ source and station from ../../../StatInfo/
StatInfo_path = '../../../StatInfo/'
station_file = StatInfo_path + 'STATION_utm_dh_4km.txt'
nRec, R, statnames = read_stat_name_utm(station_file)

source_file = StatInfo_path + 'SOURCE_SRF.txt'
nShot, S, sNames = read_source_new_utm(source_file)

# Time serie:
num_pts = 1500;
dt = 0.16;
t = np.arange(num_pts) * dt
t_min = 0;
t_max = max(t);
# Filtering parameters:
fs = 1 / dt
highcut = 0.1
lowcut = 0.025
# flo from e3d.par file for fw modeling:
flo = 1.0
delay_Time = (3 / flo)
t_shift = 0

# Name of channels
GV = ['.090', '.000', '.ver']  # emod3d output
BH = ['HHX', 'HHY', 'HHZ']  # GeoNet Broadband sensors
GV_ascii = ['.x', '.y', '.z']  # Adjoint source ascii input


GV=['.090','.000','.ver']
BH=['HHX','HHY','HHZ']

GV=['.090','.000','.ver']
GV_ascii=['.x','.y','.z']

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

# List of events included in the inversion:
#ishot_arr=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
ishot_arr=np.linspace(1,nShot,nShot).astype('int')

#Initiate the misfit sum and number of windows picked:
Err_mf_all=0
nwin_all=0
# Run data achieving, processing (shifting, tapering, filtering, integrating and formating to SAC),
# pyflex pick and pyadjoint misfit calculation:
for ishot_id in range(0,len(ishot_arr)):
    # source ishot/ event name: sNames[ishot-1]
    ishot=ishot_arr[ishot_id]
    n_ishot=0
    Err_mf_ishot=0

    Err_mf=0
    n_ishot=0
    # stations in the loop:
    for i,statname in enumerate(statnames):

        for k in range(0,3): 
            s0=statname+GV[k]
            # Read synthetic data for the updated model with the current gardients and step length:
            data_sim = timeseries.read_ascii('../../Vel_opt/Vel_ob_'+str(ishot)+'/'+s0)
            data_sim = time_shift_emod3d(data_sim,delay_Time,dt)
            
            try:
                data_obs = timeseries.read_ascii('../../Vel_ob_240s_HH_13s/Vel_ob_'+str(ishot)+'/'+s0)
            except:
                continue    # continue here if no record for this station/ component
            # Tapering:
            data_obs  = np.multiply(signal.tukey(int(num_pts),0.05),data_obs)
            data_sim  = np.multiply(signal.tukey(int(num_pts),0.05),data_sim)
            # Filtering
            data_obs = butter_bandpass_filter(data_obs, lowcut, highcut, fs, order=4)
            data_sim = butter_bandpass_filter(data_sim, lowcut, highcut, fs, order=4)
            # Integrating:
            data_obs = np.cumsum(data_obs)*dt
            data_sim = np.cumsum(data_sim)*dt
            # Formating to sac:
            obs_sac = SACTrace(kstnm=statname, kcmpnm=BH[k], stla=R[i,1], stlo=R[i,0], evla=S[ishot-1,1], evlo=S[ishot-1,0], evdp=S[ishot-1,2], nzyear=2000, nzjday=1, nzhour=0, nzmin=0, nzsec=0, nzmsec=0,
                           t0=t_min, t1=t_max, delta= dt, b=t_shift, data=data_obs)
            obs_sac.write('OBS.SAC', byteorder='little')

            syn_sac = SACTrace(kstnm=statname, kcmpnm=BH[k], stla=R[i,1], stlo=R[i,0], evla=S[ishot-1,1], evlo=S[ishot-1,0], evdp=S[ishot-1,2], nzyear=2000, nzjday=1, nzhour=0, nzmin=0, nzsec=0, nzmsec=0,
                           t0=t_min, t1=t_max, delta= dt, b=t_shift, data=data_sim)
            syn_sac.write('SYN.SAC', byteorder='little')  
            # read sac-files:
            obs_data = obspy.read("OBS.SAC")
            synth_data = obspy.read("SYN.SAC")
            # sac-processing:
            obs_data.detrend("linear")
            obs_data.detrend("demean")
            obs_data.taper(max_percentage=0.05, type="hann")
#                    obs_data.filter("bandpass", freqmin=1.0 / 20.0, freqmax=1.0 / 2.0,
#                                    corners=4, zerophase=True)
#                    
            synth_data.detrend("linear")
            synth_data.detrend("demean")
            synth_data.taper(max_percentage=0.05, type="hann")
#                    synth_data.filter("bandpass", freqmin=1.0 / 20.0, freqmax=1.0 / 2.0,
#                                      corners=4, zerophase=True)                    
            
            # Select windows
            #            windows = pyflex.select_windows(obs_data1, synth_data1, config, plot=True)
            try:
                #                windows = pyflex.select_windows(obs_data, synth_data, config, plot=False)
                windows = flexwin_new.select_windows(obs_data, synth_data, config, plot=False)
            except:
                windows = []
            if(len(windows)>0):
                if(1>0):
                    print(len(windows))
                    try:
                        #Pyadjoint for misfit calculation only
                        adj_src = pyadjoint.calculate_adjoint_source("multitaper_misfit", obs_data, synth_data,config_adj, windows,adjoint_src=False,plot=False)
                        print(adj_src.misfit)
                        Err_mf=Err_mf+adj_src.misfit
                        n_ishot = n_ishot+len(windows)
                    except:
                        print('raise Exception')                  
                
    #Everage the misfit for one event
    if(n_ishot>0):                        
        Err_mf_ishot = Err_mf/n_ishot
        nwin_all = nwin_all + n_ishot
        Err_mf_all = Err_mf_all+Err_mf_ishot
#Everage the misfit for all events
Err_FS = Err_mf_all/len(ishot_arr)
#Save the misfit for the updated model with the current gardients and step length:
f_err = open('err_opt.dat','w')
(np.float64(Err_FS)).tofile(f_err)     
print(Err_FS)
fi=open('nwin.dat','w')
(np.int64(nwin_all)).tofile(fi)
fi.close()
print(nwin_all)
    
    
