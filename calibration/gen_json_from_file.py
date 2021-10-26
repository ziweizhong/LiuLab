import numpy as np
import math
import matplotlib
import csv
import matplotlib.pyplot as plt
from random import *
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import time

#script to obtain calibration parameters for vials given a set of tests.  This then saves a file
#with the parameters for fit_eqn9 which is usable with the current evolver setup

NUM_VIALS = 16 # Number of vials we rotated through
CALIBRATION_RAW = "od_135_cal_new.txt"
MAXFEV = 100000000
PARAM_COUNT = 4

def fit_eqn9(x,A,B,C,D):
    y = A * np.log(B*x + C) + D
    return y

def fit_eqn9_inv(x, A, B, C, D):
    y = D - C * np.exp((A * (x + B)))
    return y

#Get data
real_OD = np.array([0.08,0.19,0.27,0.36,0.43,0.50,0.56,0.64,0.72,0.83,0.94,1.08,1.18,1.30,1.32,1.52])

detector_value = np.loadtxt(open(CALIBRATION_RAW, "rb"), delimiter=",", skiprows=0)

for row in range(NUM_VIALS):
    detector_value[:, row] = np.flip(np.roll(detector_value[:, row], -1-row))

param_set = np.zeros([PARAM_COUNT, NUM_VIALS])
bounds = ([-np.inf, -np.inf, -np.inf, -np.inf],
          [0, np.inf, np.inf, np.inf])

#Building subplots and getting plot limits
fig, axs = plt.subplots(4, 4, figsize=(14, 14))
ymin, ymax = (np.min(detector_value) * 0.9, 1.1 * np.max(detector_value))
xmin, xmax = (0, 1.2 * np.max(real_OD))

for vial in range(0,16):
    # Getting axes to plot in
    ax = axs[3-(vial//4), vial%4]

    #Getting X and Y values
    OD_X = real_OD
    DV_Y = detector_value[:,vial]

    # Plotting points
    ax.plot(OD_X, DV_Y, 'o')
    ax.set_title(("vial %d" % (vial)))
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    #Getting fit for this vial
    #Magic value to get the initial guess to go through both the first and last point
    offset = np.exp(-np.exp(1)) * np.ptp(DV_Y) / (1 - np.exp(-np.exp(1)))
    init = [-np.exp(1)/np.ptp(OD_X),
            -np.min(OD_X),
            offset / np.exp(-np.exp(1)),
            np.max(DV_Y) + offset]

    param, covar = curve_fit(fit_eqn9_inv, OD_X, DV_Y,
                             bounds = bounds,
                             p0 = init,
                             maxfev = MAXFEV)

    #Plotting fit
    #x_test = np.linspace(xmin, xmax, 300)
    #fit_y = fit_eqn9_inv(x_test, *param)
    #ax.plot(x_test, fit_y, '-')

    #Plotting Fit with original equation
    y_test = np.linspace(ymin, ymax, 300)
    fit_x = fit_eqn9(y_test, 1/param[0], -1/param[2], param[3]/param[2], -param[1])
    ax.plot(fit_x, y_test, '-')

    #Plotting iniital guess
    x_test = np.linspace(xmin, xmax, 300)
    init_y = fit_eqn9_inv(x_test, *init)
    ax.plot(x_test, init_y, '-')

    #Printing obtained parameter values for fit_eqn9 (not the inverse)
    param_set[:,vial] = np.array([1/param[0], -1/param[2], param[3]/param[2], -param[1]])
    print("A:% 9.2f \tB:% 9.3f \tC: % 9.3f \tD: % 9.3f" % tuple(param_set[:, vial]))


name = input("Please enter a name for this fit:\n")
f = open("od_calibration.json",'w')
f.write("{\"name\": \"%s\", \"coefficients\": ["%name)

for vial in range(NUM_VIALS):
    f.write("[%f,%f,%f,%f]"%(tuple(param_set[:, vial])))
    if vial < NUM_VIALS - 1:
        f.write(",")

current_time = time.time()

f.write('], "type": "exp", "timeFit": %f, "active": true, "params": ["od_135"]}'%current_time)

f.close()

plt.show()
