"""
A short Python script used to download GeoNet data

This requires ObsPy (obspy.org) and Python3

Information:
1) different FDSN clients:
https://docs.obspy.org/packages/obspy.clients.fdsn.html
2) filtering options:
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html?highlight=filter#obspy.core.stream.Stream.filter
3) available file formats:  
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.write.html#obspy.core.stream.Stream.write
"""
#import obspy

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from qcore import timeseries
import numpy as np
import os

import matplotlib.pyplot as plt
from scipy.signal import resample

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


def read_srf_source(source_file):

    with open(source_file, 'r') as f:
        lines = f.readlines()
    line0=lines[0].split()
    nShot=int(line0[0])
    S=np.zeros((nShot,3))
    sNames = []
    for i in range(1,nShot+1):
        line_i=lines[i].split()
        S[i-1,0]=line_i[0]
        S[i-1,1]=line_i[1]
        S[i-1,2]=line_i[2]
        sNames.append(line_i[3])

    return nShot, S, sNames


# vvv  SET PARAMETERS HERE vvv
# Station
client = "GEONET"  # GeoNet if you want NZ data
network = "NZ"  # this is usually static
#station = "BFZ"   
location = "*"  # GeoNet broadband locations are usually '10'
channel = "HH?"  # wildcard component to get 3 components
#channel = "LH?"  # wildcard component to get 3 components

# Event
# starttime = UTCDateTime("2014-01-20T02:52:45Z")  # e.g. Eketahuna M6.2
# endtime = starttime + 90  # end of waveform, X seconds are starttime
#starttime = UTCDateTime("2016-11-13T11:02:56Z")  # e.g. Kaikoura M7.8
#endtime = starttime + 500  # end of waveform, X seconds are starttime

# Processing
remove_response = True  # if you want to remove the instrument response
min_period = 10  # for bandpass filter, if 'None', no filter applied
max_period = 40

c = Client(client)
## Output
#output_format = "MSEED"  # also SAC, SEGY, SU etc.
## ^^^ SET PARAMETERS HERE ^^^
statnames_BB = ['QRZ','NNZ','WEL','DSZ','THZ','KHZ','INZ','LTZ','GVZ','OXZ']
#statnames_BB = ['NNZ']
source_file='SOURCE_SRF.txt'
nShot, S, sNames = read_srf_source(source_file)

#Number of time samples, time step and time vector defined similar to the synthetic data from emod3d:
nt=1500
dt=0.16
t = np.arange(nt)*dt

mkdir_p('Vel_ob_200s_HH_13s')
#convert m/s to cm/s
factor = 100

for ishot in range(1,nShot+1):
#for ishot in range(1,2):
    i_event = c.get_events(eventid=sNames[ishot-1]) 
    starttime_utc = str(i_event).split("\n")[1].split(' |')[0]
    
    starttime = UTCDateTime(starttime_utc)  
    endtime = starttime + nt*dt     
    
    for station in statnames_BB:
        #Create new folder for data storage:
        dir_ob = 'Vel_ob_200s_HH_13s/Vel_ob_'+str(ishot)
        mkdir_p(dir_ob)
        main_folder = dir_ob+'/'
        try:
            # Set the filename for saving the waveforms
            fid_out = "{net}_{sta}_{year}_{jday}".format(
                    net=network, sta=station, year=starttime.year, jday=starttime.julday)
            
            # Get waveform data as an Obspy Stream object
            
            st = c.get_waveforms(network=network, station=station, location=location,
                                 channel=channel, starttime=starttime, endtime=endtime,
                                 attach_response=True)
            
            # Write the raw data
            #st.write(fid_out + "_raw.{}".format(output_format), format=output_format)
            
            # Remove response ! Important, do not remove
            if remove_response:
                st.remove_response()
#            
#            # Preprocessing functionality
#            st.detrend("linear")  # detrend the data to remove any very-long period signal
#            st.detrend("demean")  # from the data
#            st.taper(max_percentage = 0.05)  # taper the ends since we cut the data
#            # Filter the data 
#            if min_period:
#                st.filter("bandpass", freqmin=1/max_period, freqmax=1/min_period)
#            st.detrend("linear")  # detrend and taper again incase filtering created any
#            st.detrend("demean")  # spurious signals
#            st.taper(max_percentage = 0.05)
            
            # print the stream object for information
            print(st)
            
            # plot the stream and save it to the current directory
            #st.plot(outfile="./{}.png".format(fid_out))
            st.plot()            
#Channel order: vx-vy-vz or 0/1/2 or E/N/Z                       
            vx = resample(st[0].data,nt)*factor
            timeseries.seis2txt(vx,dt,main_folder,station,'090')
            
            vy = resample(st[1].data,nt)*factor
            timeseries.seis2txt(vy,dt,main_folder,station,'000')

            vz = resample(st[2].data,nt)*factor
            timeseries.seis2txt(vz,dt,main_folder,station,'ver')
            plt.figure(figsize=(10,2.5))  
            plt.plot(t,vz,'r');plt.plot(t,vy,'g');plt.plot(t,vx,'b');
            plt.legend(('NZ.NNZ.10.HHE','NZ.NNZ.10.HHN','NZ.NNZ.10.HHZ'))
            plt.xlabel('Time [s]',fontsize=14)
            plt.xlim([0,240])
            plt.ylabel('cm/s',fontsize=14)
            plt.title('Event: 2011-04-29T19, Station:NNZ, Sensor: Broadband seismometer',fontsize=14)
            
            # Save the data based on the User-defined file format
            #st.write(fid_out + "_processed.{}".format(output_format), format=output_format)
        except:
            print('no record for station'+station)
    




