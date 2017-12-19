"""
- Abdelrahman Elbehery
- Section: 2
- Biomedical Engineering Assignment-2
"""
import numpy as np
import pylab
from EMG_FILTER import EMGFilter



F = EMGFilter()
# noise is chosen visually from the graph
TH = 3*np.std(np.abs(F.data[500:800]))


F.filter_avg(20, TH)
F.match_templates(threshold=TH, r_low=30000, r_high=35000)
F.plot_templates(2)