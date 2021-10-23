"""
ZZ update 2021-10-23:
This code should be used in conjunction with calibrate-stacked.py and is currently very crude.
This script will take data from outfile_vial%d and plot the first four lines, w ith the first line being x1, second line being y1, third line x2, fourth line being y2. The (x1,y1) data and (x2,y2) data will be in different colors.
This is useful in overlaying data from different calibration runs and observing the consistency of data obtained.
"""

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(4, 4)
for i in range(0,16):
    file_name = "outfile_vial%d"%i
    data = np.genfromtxt(file_name,delimiter=',')

    data1x = data[0,:15]
    data1y = data[1,:15]

    data2x = data[2,:15]
    data2y = data[3,:15]

    ax[i // 4, (i % 4)].plot(data1x, data1y, 'o', markersize=1, color='black')
    ax[i // 4, (i % 4)].plot(data2x, data2y, 'o', markersize=1, color='red')

plt.show()
