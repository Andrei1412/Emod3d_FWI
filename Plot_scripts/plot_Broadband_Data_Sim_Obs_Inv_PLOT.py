#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 15:43:17 2018
    
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy import signal
from qcore import timeseries
from scipy import integrate
from scipy.signal import butter, lfilter
#from statistics import median

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


def readGP(loc, fname):
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
    num_pts=float(line1[0])
    dt=float(line1[1])
    shift=float(line1[4])

    return data, num_pts, dt, shift

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
    num_pts=float(line1[0])
    dt=float(line1[1])
    shift=float(line1[4])

    return data, num_pts, dt, shift

def adjust_for_time_delay(ts, dt, shift):
    """
        ts: time series data
    """
    t0_index = int(shift/dt)
    if t0_index == 0:
        num_pts = ts.size
    elif t0_index > 0:
        ts = np.concatenate((np.zeros(t0_index), ts))
        num_pts = ts.size
    elif t0_index <0:
        ts = ts[np.abs(t0_index):]
        num_pts = ts.size

    return ts, num_pts, dt

def get_adjusted_stat_data(loc, stat_code):
    """
    Time series data is adjust for time shift

    returns a dictionary with components:
        000, 090, ver, t, stat_code
    """
    stat_data = {"000": None, "090": None, "ver":None, "t": None, "name": None}
    g=981. #cm/s^2
    stat_data["000"], num_pts, dt, shift  = readGP(loc, ".".join([stat_code, "000"]))
    stat_data["090"], num_pts, dt, shift = readGP(loc, ".".join([stat_code, "090"]))   
    stat_data["ver"], num_pts, dt, shift = readGP(loc, ".".join([stat_code, "ver"]))


    stat_data["000"], num_pts, dt = adjust_for_time_delay(stat_data["000"], dt, shift)
    stat_data["090"], num_pts, dt = adjust_for_time_delay(stat_data["090"], dt, shift)
    stat_data["ver"], num_pts, dt = adjust_for_time_delay(stat_data["ver"], dt, shift)

    t = np.arange(num_pts)*dt
    stat_data["t"] = t

    stat_data["name"]=stat_code
    return stat_data

def read_stat_name(station_file):

    with open(station_file, 'r') as f:
        lines = f.readlines()
    line0=lines[0].split()
    nRec=int(line0[0])
    R=np.zeros((nRec,3))
    statnames = [] 
    for i in range(1,nRec+1):
        line_i=lines[i].split()
        R[i-1,0]=int(line_i[0])
        R[i-1,1]=int(line_i[1])
        R[i-1,2]=int(line_i[2])
        statnames.append(line_i[3])
    return nRec, R, statnames

def rms(stat_data):

    num_pts=len(stat_data)
    D = (np.sum(np.square(stat_data))/num_pts)**0.5

    return stat_data/D


def ncc(stat_data_0_Sf,stat_data_0_Of,num_pts, delta_T,dt):
    """
    Normalized correlation coefficient
    """    
#    rwm1=0
#    rwm2=0
#    rwm3=0
    num_delta_t = int(delta_T/dt);    td = np.arange(-num_delta_t,num_delta_t)*dt;
    t = np.arange(num_pts)*dt
    ncc_array=np.zeros(len(td))
    
    for it in range(0,len(td)):
        stat_data_0_S_shift = np.zeros(num_pts); n_shift=int(np.abs((td[it]/dt)));
        
        if td[it]<0:
            stat_data_0_S_shift[n_shift:num_pts] = stat_data_0_Sf[0:num_pts-n_shift]                             
        else:
            stat_data_0_S_shift[0:num_pts-n_shift] = stat_data_0_Sf[n_shift:num_pts]                     
        
        rwm1_arr=np.multiply(stat_data_0_S_shift,stat_data_0_Of)
        rwm2_arr=np.square((stat_data_0_S_shift))
        rwm3_arr=np.square((stat_data_0_Of))        
        
        rwm1=integrate.simps(rwm1_arr,t)
        rwm2=integrate.simps(rwm2_arr,t)
        rwm3=integrate.simps(rwm3_arr,t)
    
        ncc_array[it]=rwm1/((rwm2*rwm3)**0.5)
    
    ncc_max = np.max(ncc_array); id_max = np.argmax(ncc_array); td_max = td[id_max];
    
    return ncc_max, td_max  

def time_shift_emod3d(data,delay_Time,dt):
    n_pts = len(data)
    ndelay_Time = int(delay_Time/(dt))
    data_shift = np.zeros(data.shape)
    data_shift[0:n_pts-ndelay_Time] = data[ndelay_Time:n_pts]
    return data_shift 

station_file = 'STATION_dh_4km.txt'  

nRec, R, statnames = read_stat_name(station_file)
statnames=statnames[0:4]
print('statnames')
print(statnames)

num_pts = 1500;
dt = 0.16;
t = np.arange(num_pts) * dt

############/nesi/nobackup/nesi00213/RunFolder/tdn27/rgraves/Adjoint/Syn_VMs/Kernels/#########################
fs = 1/dt
lowcut = 0.05
highcut = 0.1

fc = highcut  # Cut-off frequency of the filter
w = fc / (fs / 2) # Normalize the frequency
b, a = signal.butter(4, w, 'low')

GV=['.090','.000','.ver']
vxyz=['vx','vy','vz']

ishot = 1
os.system('/scale_wlg_nobackup/filesets/nobackup/nesi00213/RunFolder/tdn27/rgraves/NZVMs/Marlborough_Events_4KM/NewVM_20200207/INV_MarlVM_DH_4km_template/Kernels/Vel_ob_240s_HH_13s/Vel_ob_'+str(ishot)+'/*.* /home/andrei/workspace/GMPlots/Sim/Vel_ob/')     
os.system('/scale_wlg_nobackup/filesets/nobackup/nesi00213/RunFolder/tdn27/rgraves/NZVMs/Marlborough_Events_4KM/NewVM_20200207/INV_MarlVM_DH_4km_template/Kernels/Vel_es/Vel_es_'+str(ishot)+'/*.* /home/andrei/workspace/GMPlots/Sim/Vel_es/')     
os.system('/scale_wlg_nobackup/filesets/nobackup/nesi00213/RunFolder/tdn27/rgraves/NZVMs/Marlborough_Events_4KM/NewVM_20200207/INV_MarlVM_DH_4km_template/Kernels/Vel_opt/Vel_ob_'+str(ishot)+'/*.* /home/andrei/workspace/GMPlots/Sim/Vel_es_inv/')
# flo from e3d.par file for fw modeling:
flo = 1.0
delay_Time = (3 / flo)
t_shift = 0


ncc_R = np.zeros([nRec,3]) 
td_max_R = np.zeros([nRec,3]) 

delta_T = 5 #5s

for i,statname in enumerate(statnames):    
        for k in range(0,3): 
            s0=statname+GV[k]
            # Read synthetic data for the current model at i-th itearation:
            data_sim = timeseries.read_ascii('/home/andrei/workspace/GMPlots/Sim/Vel_es/' + s0)
            data_sim = time_shift_emod3d(data_sim, delay_Time, dt)
            
            data_sim_inv = timeseries.read_ascii('/home/andrei/workspace/GMPlots/Sim/Vel_es_inv/' + s0)
#            data_sim_inv = timeseries.read_ascii('/home/andrei/workspace/GMPlots/Sim/Vel_es/' + s0)            
            data_sim_inv = time_shift_emod3d(data_sim_inv, delay_Time, dt)            

            try:
                data_obs = timeseries.read_ascii('/home/andrei/workspace/GMPlots/Sim/Vel_ob/' + s0)
            except:
                continue  # continue here if no record for this station/ component
            # Tapering:
            data_obs = np.multiply(signal.tukey(int(num_pts), 0.05), data_obs)
            data_sim = np.multiply(signal.tukey(int(num_pts), 0.05), data_sim)
            data_sim_inv = np.multiply(signal.tukey(int(num_pts), 0.05), data_sim_inv)            
            # Filtering
            data_obs = butter_bandpass_filter(data_obs, lowcut, highcut, fs, order=4)
            data_sim = butter_bandpass_filter(data_sim, lowcut, highcut, fs, order=4)          
            data_sim_inv = butter_bandpass_filter(data_sim_inv, lowcut, highcut, fs, order=4)               
            
            fft_obs=fftpack.fft(data_obs)
            fft_sim=fftpack.fft(data_sim)
            fft_sim_inv=fftpack.fft(data_sim_inv)
            
            fft_obs=np.abs(fft_obs)
            fft_sim=np.abs(fft_sim)
            fft_sim_inv=np.abs(fft_sim_inv)
            
            sample_freq = fftpack.fftfreq(data_obs.size, d=dt)
            
            plt.figure(figsize=(10,2.5))       
            plt.subplot(1,2,1)
        #    plt.figure(figsize=(10,2.5))    
            plt.plot(t,data_obs,c='k',linestyle='solid')
            plt.plot(t,data_sim,c='r',linestyle='dashed')
            plt.plot(t,data_sim_inv,c='r',linestyle='solid')
        
            plt.title('Observed vs Simulated seismograms at station '+statname, loc='center')
            plt.xlabel('Time (s)')
            plt.ylabel(vxyz[k]+' (cm/s)')
            plt.gca().legend(('Observed','Init.Simulated','Inv.Simulated'))
        #    plt.xlim([0,100])
            #ylimit = 1.1*max([ymax1,ymax4])
        #    plt.plot([0,20],[-ylimit*0.50,-ylimit*0.50],c='k',linewidth=1.0)
        #    plt.vlines([0,20],-ylimit*0.6,-ylimit*0.4,color='k',linewidth=1.0)
#            plt.ylim([-ylimit,ylimit])
            ax = plt.gca()
            #ax.text(0.05,0.2, '20s', horizontalalignment='left', verticalalignment='top', transform = ax.transAxes, fontsize=9.0)
            ax.axis('on')
            plt.grid()
            plt.subplot(1,2,2)
        #    plt.figure(figsize=(10,2.5))    
            plt.plot(sample_freq,fft_obs,c='k',linestyle='solid')
            plt.plot(sample_freq,fft_sim,c='r',linestyle='dashed')
            plt.plot(sample_freq,fft_sim_inv,c='r',linestyle='solid')
            plt.title('FFT', loc='center') 
            plt.xlim([0,.2])
#            plt.ylim([0,yflimit])    
            plt.grid()
            plt.xlabel('Frequency (Hz)')
            plt.gca().legend(('Observed','Init.Simulated','Inv.Simulated'))
            plt.subplots_adjust(left=.1, bottom=0.0, right=2.0, top=1.0, wspace=0.2, hspace=-0.3)
            #plt.savefig(statname + '.png',dpi=200)
            plt.show()
        
            ncc_R[i,k], td_max_R[i,k] = ncc(fft_sim,data_obs,num_pts, delta_T,dt)
            
