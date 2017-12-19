"""
Name : Abdelrahman Elbehery
Section : 2
ID : 1300759
Assignment : 1- ECG detector
"""

import numpy as np
import pylab

class EMGFilter(object):
    """
    A class implementing the EMG filter. The filter performs the following sequence of operations
    + Rectifying the signal by removing -ve samples [absolute operation]
    + Then an average window algorithm is applied to the samples for a smoother curve
    + A threshold is assigned to capture the peaks to detect MUAPs
    + The matching algorithm is invoked to capture the repeated templates of MUAPs

    ## Attributes
    + data : np.array that containts the sampled ECG data,
    exist only for visualization and debugging
    + data_filtered_avg : np.array that contains the sampled filtered data after applying
    moving window with size N
    + r_peaks : np.array that contains the peak indeces of the EMG
    + r_peaks_idx : np.array that contains the index number of each peak
    + template_matrix : numpy 2d array where each row represent a detected template
    a vectorized difference operation is done over this array to match other templates
    + template_count : a numpy 1D array that holds the number of occurrences of each template.
    This attribute is used to check if such a template is being repeated often or
    it is a superposition. Also a zero value means we dont have to draw such patteren
    + threshold : an autmatically generated value used as the signal threshold to detect peaks


    ### The features implemented are
    + Moving average technique with N as an argument
    + R-peak detection algorithm

    ## Args
    + file_name : the file name containing the ECG samples
    """

    def __init__(self, file_name="Data.txt"):
        self.data = np.loadtxt(file_name)
        self.data_filtered_avg = None
        self._data_filtered = None
        self.threshold = None
        self.r_peaks = None
        self.r_peaks_idx = None
        self.r_peaks_clrs = None
        self.template_matrix = None
        self.template_count = None

    def filter_avg(self, avg_window_size, threshold):
        """
        Attempts to rectify the signal and perform the vectorized moving average algorithm
        followed by the peak detection algorithm.
        The method exists only for debugging and visualization, the user should use
        the get_templates() interface that will automatically invoke this method
        to do the proper pre-processing for the template matching algorithm

        The method mutates the following attributes
        + _data_filtered : an internal variable that has the rectified data
        + data_filtered_avg : the data after the moving average algorithm
        + r_peaks : an array of non-zero values at the peak indices
        + r_peaks_idx : an array of the indices of the peaks detected

        ## Args
        + avg_window_size : the average window number of samples
        + threshold : threshold for peak detection algorithm @ the moving average.
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
            # only compare it to the past value, as long as it exists distance > T
            val = val * (val > self.threshold) *\
            (val > past_val or (max_idx - past_max_idx) > avg_window_size)
            past_val = past_val * (past_val > val)
            self.r_peaks[max_idx] = val
            self.r_peaks[past_max_idx] = past_val
        self.r_peaks_idx = np.nonzero(self.r_peaks)[0]#.astype('float')
        # dont show zeros, useful for plotting
        self.r_peaks[self.r_peaks == 0] = np.nan


    def match_templates(self, avg_window_size=20, threshold=11.7, diff_th=12.65**5):
        """
        Capture the different templates of the MUAPs, in case a template is repeated for one or two
        times only the algorithm considers it as a superposition happened between multiple MUAPs.
        It will still be shown, but is labled by "superposition"

        ## Algorithm
        The function relies heavily on the vectorized matrix operations provided by numpy.
        Since it does the difference across all the matrix elements by broadcasting the last
        captured template. The resulted value is a vector of distances, the argmax of this distance
        is considered the winning template. Once this winning template is known and is updated
        to the last captured template

        ## Arguments
        + avg_window_size : the average window number of samples
        + threshold : threshold for peak detection algorithm @ the moving average. **Note**
        at this point, the theshold is calculated manually by inspecting the signal
        + diff_th : difference threshold for template matching
        """
        # filter the data, obtain the peak indices
        self.filter_avg(avg_window_size, threshold)
        # now create the template matrix with a total number of templates = total number of peaks
        # this might not hold true, since many of these peaks are repeated. But we are allocating
        # an upper bound size for this numpy array
        self.template_matrix = np.zeros([self.r_peaks_idx.shape[0], avg_window_size])
        self.template_count = np.zeros([self.r_peaks_idx.shape[0]]).astype('int32')
        self.r_peaks_clrs = np.zeros([self.r_peaks_idx.shape[0], 3])
        next_temp_pos = 0
        # now loop through all the detected peaks, and compate it [vector difference] against
        # all the detected templates
        for p_idx, t_idx in enumerate(self.r_peaks_idx):
            dist = int(avg_window_size//2)
            template = self.data_filtered_avg[t_idx-dist:t_idx+dist]
            #dist_vect = np.linalg.norm(self.template_matrix - template, axis=1) < diff_th
            dist_vect = np.sum((self.template_matrix - template)**2, axis=1) < diff_th
            win_idx = np.argmax(dist_vect)
            # if t_idx==30299 or t_idx==30688:
            #     pass
            if dist_vect[win_idx] == 0 or self.template_count[win_idx] == 0:
                # None of the templates matched this one, add it into the template matrix
                self.template_matrix[next_temp_pos, :] = template
                # Now declare that we have found one instance of such template
                self.template_count[next_temp_pos] = 1
                # assign some random color to these group of points
                clr = np.random.random(3)
                self.r_peaks_clrs[win_idx, :] = clr
                self.r_peaks_clrs[p_idx, :] = clr
                # increment the next position pointer
                next_temp_pos += 1
                continue
            # we have found a winner, increment the number of occurrences
            self.template_count[win_idx] += 1
            self.r_peaks_clrs[p_idx, :] = self.r_peaks_clrs[win_idx, :]
            # then update this template to be the new one [SLIDE 13]
            self.template_matrix[win_idx, :] = template

    def plot_templates(self, num_fig_h=4):
        """
        Plots all the detected patterns using pylab

        ## Args
        + num_fig_h : the number of figures to draw horizontaly
        """
        idxs = np.nonzero(self.template_count > 0)[0]
        sub_plot_x = idxs.shape[0]//num_fig_h
        sub_plot_y = num_fig_h

        for i in range(len(idxs)):
            pylab.subplot(sub_plot_x, sub_plot_y, i+1)
            pylab.plot(self.template_matrix[idxs[i]])
            pylab.title("Template repeated: "+ str(self.template_count[idxs[i]])+ " times", loc='left')
        pylab.show()

    def plot_peaks(self, r_low=0, r_high=0):
        """
        Plots all the detected peaks while maintaining the same color
        for the peaks that belong to the same template
        """
        if (r_low != 0 and r_high != 0):
            pylab.plot(self.data_filtered_avg[r_low:r_high])
        else:
            pylab.plot(self.data_filtered_avg)

        for i, p_idx in enumerate(self.r_peaks_idx):
            if (r_low != 0 and r_high != 0) and (p_idx > r_high or p_idx < r_low):
                continue
            pylab.plot(p_idx-r_low, self.r_peaks[p_idx], '*',c=self.r_peaks_clrs[i])
        pylab.show()