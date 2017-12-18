"""
Name : Abdelrahman Elbehery
Section : 2
ID : 1300759
Assignment : 1- ECG detector
"""

import numpy as np
from scipy import signal


class EMG_Filter(object):
    """
    A class implementing the EMG filter. The filter performs
    + Notch filter to remove the power-line noise followed by a BPF
    + A 5 point difference to check large slope regions
    + Then square the difference array to poroduce positive samples
    + Then an average window algorithm is applied to the samples for a smoother curve
    + A threshold is assigned to capture the QRS [the R peak]

    ## Attributes
    + data : np.array that containts the sampled ECG data, exist only for visualization
    + data_filtered_avg : np.array that contains the sampled filtered data after applying
    moving window with size N
    + r_peaks : np.array that contains the peak indeces of the ECG [R portion]
    + time_stamp : np.array that contains the time of each R peak
    + threshold : an autmatically generated value used as the signal threshold to detect peaks


    ### The features implemented are
    + Moving average technique with N as an argument
    + Threshold value calculation
    + R-peak detection algorithm

    ## Args
    + file_name : the file name containing the ECG samples
    + f_sampling : the sampling frequency
    """

    def __init__(self, file_name="Data.txt"):
        self.data = np.loadtxt(file_name)
        self.data_filtered_avg = None
        self._data_filtered = None
        self.r_peaks = None
        self.threshold = None

    def filter_avg(self, avg_window_size, threshold):
        """
        Performs all the filter function on the loaded data. The method will not
        perform the first two filtering stages [finite diff, square] if dont_filter==True

        ## Args
        + avg_window_size : the average window number of samples
        """

        # perform absolute rectification
        self._data_filtered = np.abs(self.data)

        # moving average
        avg_window_size -= 1
        self.data_filtered_avg = np.copy(self._data_filtered[avg_window_size:])
        for i in range(avg_window_size):
            self.data_filtered_avg += self._data_filtered[i:-avg_window_size+i]
        self.data_filtered_avg /= avg_window_size

        self.r_peaks = np.zeros([self.data_filtered_avg.shape[0],])
        self.threshold = threshold
        ## peak detection
        for i in range(0, self.data_filtered_avg.shape[0], avg_window_size):
            max_idx = np.argmax(self.data_filtered_avg[i:i+avg_window_size])+i
            if i >= avg_window_size:
                past_max_idx = np.argmax(self.r_peaks[i-avg_window_size:i])+(i-avg_window_size)
            else:
                past_max_idx = 0
            val = self.data_filtered_avg[max_idx]
            past_val = self.r_peaks[past_max_idx]
            val = val * (val > self.threshold) * (val > past_val)
            past_val = past_val * (past_val > val)
            self.r_peaks[max_idx] = val
            self.r_peaks[past_max_idx] = past_val
        #self.rr_interval = np.nonzero(self.r_peaks)[0].astype('float')
        #self.time_stamp = self.rr_interval#*1/self.f_sampling
        # dont show zeros
        self.r_peaks[self.r_peaks == 0] = np.nan
        