"""
Name : Abdelrahman Elbehery
Section : 2
ID : 1300759
Assignment : 1- ECG detector
"""

import numpy as np
from scipy import signal
import pylab


class Filter(object):
    """
    A class impleminting the ECG filter. The filer performs
    + A 5 point different to check large slope regions
    + Then square the difference array to poroduce positive samples
    + Then an average window algorithm is applied to the samples for a smoother curve
    + A threshold is assigned to capture the QRS [the R peak]

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
        ### A built-in function implements the filter
        w0 = notch_freq/(self.f_sampling/2)
        b, a = signal.iirnotch(w0, 30.0)
        self.data_filtered = signal.lfilter(b, a, self.data)
        #####################################################

        self.data_filtered = (self.f_sampling/(8.0)*(-self.data_filtered[2:-2]-\
        2*self.data_filtered[1:-3] + 2*self.data_filtered[3:-1] + self.data_filtered[4:]))**2

        avg_window_size -= 1
        self.data_filtered_avg = np.copy(self.data_filtered[avg_window_size:])
        for i in range(avg_window_size):
            self.data_filtered_avg += self.data_filtered[i:-avg_window_size+i]

        self.data_filtered_avg /= avg_window_size
   


x=Filter()
x.filter_avg(20)
pylab.subplot(2,1,1)
pylab.plot(x.data[:3000])
pylab.subplot(2,1,2)
pylab.plot(x.data_filtered_avg[:3000])
pylab.show()
