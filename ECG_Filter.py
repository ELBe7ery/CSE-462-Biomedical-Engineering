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
    + Threshold detection

    ## Args
    + file_name : the file name containing the ECG samples
    + threshold : the threshold needed to detect the R peak. If not given it will be calculated
    by the following formula (max(ecg)-avg(ecg)/2)
    + f_sampling : the sampling frequency
    """

    def __init__(self, file_name="DataN.txt", threshold=None, f_sampling=256.0):
        self.data = np.loadtxt(file_name)
        # for display purposes
        self.data_filtered = None
        # for display purposes
        self.data_filtered_avg = None
        self._data_filtered = None
        self.f_sampling = f_sampling
        self.threshold = threshold

    def filter_avg(self, avg_window_size, notch_freq=50.0):
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

