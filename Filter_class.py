"""
Name : Abdelrahman Elbehery
Section : 2
ID : 1300759
Assignment : 1- ECG detector
"""

import numpy as np
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

    def filter_avg(self, avg_window_size):
        """
        Performs all the filter function on the loaded data. The method will not
        perform the first two filtering stages [finite diff, square] if dont_filter==True

        ## Args
        + avg_window_size : the average window number of samples
        """
        self.data_filtered = (self.f_sampling/(8.0)*(-self.data[2:-2]- 2*self.data[1:-3]\
        + 2*self.data[3:-1] + self.data[4:]))**2

        avg_window_size-=1
        self.data_filtered_avg = np.copy(self.data_filtered[avg_window_size:])
        for i in range(avg_window_size):
            self.data_filtered_avg += self.data_filtered[i:-avg_window_size+i]

        self.data_filtered_avg/=avg_window_size
    


x=Filter()
x.filter_avg(50)
pylab.subplot(2,1,1)
pylab.plot(x.data_filtered[:1500])
pylab.subplot(2,1,2)
pylab.plot(x.data_filtered_avg[:1500])
pylab.show()

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

DATA = np.loadtxt("DataN.txt")
#DATA = running_mean(DATA,50)
#pylab.plot(DATA[2000:])
#pylab.show()
# #pylab.plot(np.diff(DATA, 256))

# pylab.subplot(2,1,1)
# pylab.plot(DATA[:2000])

# pylab.subplot(2,1,2)
    # DATA = 256/(8.0)*(-DATA[2:-2] - 2*DATA[1:-3] + 2*DATA[3:-1] + DATA[4:])
    # DATA = DATA**2

    # pylab.subplot(2,1,2)
    # DATA = running_mean(DATA,50)
    # pylab.plot(DATA[:2000])
    # pylab.show()

# pylab.plot((DATA**2)[:2000])

# pylab.show()