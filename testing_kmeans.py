"""
Testing k means
"""
from sklearn.cluster import KMeans
from EMG_FILTER import EMGFilter
import numpy as np
import pylab
np.random.seed(12)
UUT = EMGFilter()
# noise is chosen visually from the graph
TH = 2000 #3*np.std(np.abs(UUT.data[500:800]))
UUT.match_templates_kmeans(threshold=TH, num_clusters=2, low=30000, high=35000)
UUT.plot_templates(2, 1, r_low=30000, r_high=35000)
UUT.plot_peaks(r_low=30000, r_high=35000)
#pylab.plot(UUT.data_filtered_avg)
#pylab.show()
# UUT.filter_avg(20, TH)
# DATA_SET = np.zeros([UUT.r_peaks_idx.shape[0], 20])

# for p_idx, t_idx in enumerate(UUT.r_peaks_idx):
#     #print (t_idx)
#     DATA_SET[p_idx, :] = UUT.data_filtered_avg[t_idx-10:t_idx+10]
# kmeans = KMeans(n_clusters=2, random_state=0).fit(DATA_SET)
# #kmeans.predict(DATA_SET[0])