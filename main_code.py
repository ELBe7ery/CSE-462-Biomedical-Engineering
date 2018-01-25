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
    UUT.match_templates(threshold=TH, low=30000, high=35000)
    UUT.plot_templates(1)
    UUT.plot_peaks(r_low=30000, r_high=35000)
question_1()