"""
- Abdelrahman Elbehery
- Section: 2
- Biomedical Engineering Assignment-2
"""
import numpy as np
import pylab
from EMG_FILTER import EMGFilter


def question_1():
    """
    Solves Q-1
    """
    UUT = EMGFilter()
    # noise is chosen visually from the graph
    TH = 3*np.std(np.abs(UUT.data[500:800]))
    UUT.filter_avg(20, TH)
    UUT.match_templates(threshold=TH)
    UUT.plot_templates(2)
    UUT.plot_peaks(r_low=30000, r_high=35000)


question_1()