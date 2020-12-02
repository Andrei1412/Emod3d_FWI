import numpy as np
from scipy import signal
from scipy.signal import butter, lfilter
from qcore import timeseries

def winpad(lt, t_off, t_on, pad):
    #window by lt:length of the time serie, t_off:end time (index), t_on: start time (index)
    """
    Sin taper window
    """
    # pad=5
    L = t_off - t_on + 2 * pad
    window = np.ones((L))

    # x=np.arange(0,pad,1)
    x = np.linspace(0, np.pi / 2, pad)
    sinx = np.sin(x)
    window[0:pad] = sinx
    window[L - pad:L] = sinx[::-1]
    print('lt=' + str(lt))
    ar1 = np.zeros((t_on - pad))
    ar2 = np.zeros((lt - t_off - pad))
    window_pad0 = np.concatenate((ar1, window))
    window_pad = np.concatenate((window_pad0, ar2))

    return window_pad


############################################################

def butter_bandpass(lowcut, highcut, fs, order):
    #Butter bandpass filter coeeficients from scipy.signal
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    #Butter bandpass filter from scipy.signal with data
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def rms(stat_data):
    #relative waveform misfit
    num_pts = len(stat_data)
    D = (np.sum(np.square(stat_data)) / num_pts) ** 0.5

    return stat_data / D


def time_shift(data, t_shift, dt):
    #shift data forward
    n_pts = len(data)
    nshift_T = int(t_shift / (dt))
    data_shift = np.zeros(data.shape)
    data_shift[nshift_T:n_pts] = data[0:n_pts - nshift_T]
    return data_shift


def time_shift_emod3d(data, delay_Time, dt):
    #shift data backward (to correct for emod3d flo-padding of delay_Time = 3/flo)
    n_pts = len(data)
    ndelay_Time = int(delay_Time / (dt))
    data_shift = np.zeros(data.shape)
    data_shift[0:n_pts - ndelay_Time] = data[ndelay_Time:n_pts]
    return data_shift


def write_win_qual(R_Time_record_end, dt, windows, filename):
    #write window from pyflex in flewin format
    fid = open(filename, 'w')
    fid.write("%s\n" % ('# NUM_WIN ='))
    fid.write("%s\n" % ('# i win_start win_end Tshift CC dlnA'))

    for i in range(0, len(windows)):
        t_off = windows[i].right * dt
        if (R_Time_record_end < windows[i].right * dt):
            t_off = R_Time_record_end
        fid.write("%4d%20f%20f%20f%20f%20f\n" % (
        i + 1, windows[i].left * dt, t_off, windows[i].cc_shift * dt, windows[i].max_cc_value, windows[i].dlnA))
    return


def write_adj_source_ts(v1, mainfolder_source, source, dt):
    #write adjoint source/ ascii file of name A.x e.g.
    # filename1=mainfolder_source+v1
    vs1 = v1.split('.')
    timeseries.seis2txt(source, dt, mainfolder_source, vs1[0], vs1[1])
    return
