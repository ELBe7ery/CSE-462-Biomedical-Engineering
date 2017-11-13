"""
Name : Abdelrahman Elbehery
Section : 2
ID : 1300759
Assignment : 1- ECG detector
"""

import numpy as np
from scipy import signal
import pylab


class ECG_Filter(object):
    """
    A class impleminting the ECG filter. The filter performs
    + Notch filter to remove the power-line noise followed by a BPF
    + A 5 point difference to check large slope regions
    + Then square the difference array to poroduce positive samples
    + Then an average window algorithm is applied to the samples for a smoother curve
    + A threshold is assigned to capture the QRS [the R peak]

    ### The only used built-in functions are
    + sympy.signal.iirnotch : to generate the notch filter components
    + sympy.signal.butter : to generate the band pass filter components
    + sympy.signal.lfilter : to apply the filter over the signals

    ### The features implemented are
    + 5 Point difference algorithm
    + Moving average technique with N as an argument
    + Threshold value calculation
    + R-peak detection algorithm

    ## Args
    + file_name : the file name containing the ECG samples
    + threshold : the threshold needed to detect the R peak. If not given it will be calculated
    by the following formula (max(ecg)-avg(ecg)/2)
    + f_sampling : the sampling frequency
    """

    def __init__(self, file_name="data_set//DataN.txt", f_sampling=256.0):
        self.data = np.loadtxt(file_name)
        # for display purposes
        self.data_filtered = None
        # for display purposes
        self.data_filtered_avg = None
        self._data_filtered = None
        self.r_peaks = None
        self.rr_interval = None
        self.f_sampling = f_sampling
        self.threshold = None

    def filter_avg(self, avg_window_size, notch_freq=50.0, no_filter=0):
        """
        Performs all the filter function on the loaded data. The method will not
        perform the first two filtering stages [finite diff, square] if dont_filter==True

        ## Args
        + avg_window_size : the average window number of samples
        + notch_freq : desired freq. to block by the notch filter
        """
        #####################################################
        ### NOTCH FILTER REGION
        ### A built-in function that implements the filter
        w0 = notch_freq/(self.f_sampling/2)
        b, a = signal.iirnotch(w0, 30.0)
        self.data_filtered = signal.lfilter(b, a, self.data)
        nyq = 0.5 * self.f_sampling
        low = 0.1 / nyq
        high = 45 / nyq
        b, a = signal.butter(2, [low, high], btype='band')
        self.data_filtered = signal.lfilter(b, a, self.data_filtered)
        if no_filter==1 : self.data_filtered = self.data
        #####################################################

        # perform 5-point difference
        self._data_filtered = (self.f_sampling/(8.0)*(-self.data_filtered[2:-2]-\
        2*self.data_filtered[1:-3] + 2*self.data_filtered[3:-1] + self.data_filtered[4:]))**2

        # moving average
        avg_window_size -= 1
        self.data_filtered_avg = np.copy(self._data_filtered[avg_window_size:])
        for i in range(avg_window_size):
            self.data_filtered_avg += self._data_filtered[i:-avg_window_size+i]
        self.data_filtered_avg /= avg_window_size

        self.r_peaks = np.zeros([self.data_filtered_avg.shape[0],])
        th_ratio = 0.6
        if no_filter==1:
            # bigger ratio to avoid capturing low points
            th_ratio = 0.8
        self.threshold = th_ratio*np.max(self.data_filtered_avg)
        for i in range(0, self.data_filtered_avg.shape[0], avg_window_size):
            max_idx = np.argmax(self.data_filtered_avg[i:i+avg_window_size])+i
            if i >= avg_window_size:
                past_max_idx = np.argmax(self.r_peaks[i-avg_window_size:i])+(i-avg_window_size)
            else:
                past_max_idx=0
            val = self.data_filtered_avg[max_idx]
            past_val = self.r_peaks[past_max_idx]
            val = val * (val > self.threshold) * (val > past_val)
            past_val = past_val * (past_val > val)
            self.r_peaks[max_idx] = val
            self.r_peaks[past_max_idx] = past_val
        self.rr_interval = np.nonzero(self.r_peaks)[0].astype('float')
        # dont show zeros
        self.r_peaks[self.r_peaks == 0] = np.nan
        self.rr_interval = self.rr_interval[1:] - self.rr_interval[:-1]
        # msec scale
        self.rr_interval *= 1/self.f_sampling

