"""
Name : Abdelrahman Elbehery
Section : 2
ID : 1300759
Assignment : 1- ECG detector - Question [1]
"""

from ECG_Filter import ECG_Filter
import pylab


def question_1(point):
    """
    Solves Q-1 point by point by drawing the proper figures
    """
    UUT = ECG_Filter()
    if point == 1:
        UUT.filter_avg(0)
        # plot before filter
        pylab.subplot(2, 1, 1)
        pylab.plot(UUT.data[:2000])
        pylab.xlabel("Sample")
        pylab.ylabel("Amplitude")
        pylab.title("Before Filter", loc='left')
        # plot After filter
        pylab.subplot(2, 1, 2)
        pylab.plot(UUT.data_filtered[:2000], color='red')
        pylab.xlabel("Sample")
        pylab.ylabel("Amplitude")
        pylab.title("After Filter [Notch filter + BPF]", loc='left')
        pylab.show()

    elif point == 2:
        n = [5, 15, 25]
        for i in n:
            UUT.filter_avg(i)
            pylab.plot(UUT.data_filtered_avg[:2000])
            pylab.plot(UUT.r_peaks[:2000], 'b*', markersize=7)
            pylab.title("R waves with N = "+str(i))
            pylab.show()
    elif point == 6:
        UUT.filter_avg(25, no_filter=1)
        pylab.plot(UUT.data_filtered_avg[:2000])
        pylab.plot(UUT.r_peaks[:2000], 'b*', markersize=7)
        pylab.title("NON FILTERED R waves with N = 25")
        pylab.show()
    elif point == 7:
        UUT.filter_avg(25)
        pylab.plot(UUT.rr_interval*(1e3))   # ms scale
        pylab.title("R:R interval graph")
        pylab.xlabel("Beat number")
        pylab.ylabel("RR intervals [ms]")
        pylab.ylim(0, 1170)
        pylab.xlim(0, 270)
        pylab.show()


#question_1(7)
UUT = ECG_Filter('data_set/Data2.txt')
# UUT.filter_avg(0)
# pylab.subplot(2, 1, 1)
# pylab.plot(UUT.data[:2000])
# pylab.xlabel("Sample")
# pylab.ylabel("Amplitude")
# pylab.title("Before Filter", loc='left')
# # plot After filter
# pylab.subplot(2, 1, 2)
# pylab.plot(UUT.data_filtered[:2000], color='red')
# pylab.xlabel("Sample")
# pylab.ylabel("Amplitude")
# pylab.title("After Filter [Notch filter + BPF]", loc='left')
# pylab.show()


UUT = ECG_Filter('data_set/Data2.txt')
i=25
UUT.filter_avg(i)
UUT.detect_sinusArrest()