"""
Name : Abdelrahman Elbehery
Section : 2
ID : 1300759
Assignment : 1- ECG detector - Question [1]
"""

from ECG_Filter import ECG_Filter
import pylab

UUT = ECG_Filter()
UUT.filter_avg(10)
## Q-1-a
# plot before filter
pylab.subplot(4, 1, 1)
pylab.plot(UUT.data[:2000], label="Before Filter")
pylab.legend(loc='upper right')
# plot After filter
pylab.subplot(4, 1, 2)
pylab.plot(UUT.data_filtered[:2000], label="After Filter [Notch filter + BPF]", color='red')
pylab.legend(loc='upper right')
# plot after derivative squared
# pylab.subplot(4, 1, 3)
# pylab.plot(UUT._data_filtered[:2000], label="After Filter [5-point diff squared]", color='purple')
# pylab.legend(loc='upper right')

UUT.filter_avg(25)
pylab.subplot(4, 1, 3)
pylab.plot(UUT.r_peaks[:2000], 'b*', markersize=7)
pylab.plot(UUT.data_filtered_avg[:2000], label="All filters\nWith R peaks N=25", color='purple')
pylab.legend(loc='upper right')

## Q-1-b

UUT.filter_avg(50)
pylab.subplot(4, 1, 4)
pylab.plot(UUT.r_peaks[:2000], 'b*', markersize=7)
#pylab.subplot(4, 1, 4)
pylab.plot(UUT.data_filtered_avg[:2000], label="All filters\nWith R peaks N=50", color='green')
pylab.legend(loc='upper right')
pylab.show()