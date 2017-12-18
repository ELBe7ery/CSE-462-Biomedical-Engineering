"""
- Abdelrahman Elbehery
- Section: 2
- Biomedical Engineering Assignment-2
"""
import numpy as np
import pylab


from EMG_FILTER import EMG_Filter

MY_FILTER = EMG_Filter()
THRESHOLD = 3*np.std(np.abs(MY_FILTER.data[500:800]))

MY_FILTER.filter_avg(20, THRESHOLD)
pylab.plot(MY_FILTER.data_filtered_avg[30000:35000])
pylab.plot(MY_FILTER.r_peaks[30000:35000], '*')
pylab.show()