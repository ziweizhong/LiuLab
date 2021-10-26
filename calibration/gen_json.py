import numpy as np
import math
import matplotlib
import csv
import matplotlib.pyplot as plt
from random import *
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import time

def fit_eqn9(x,A,B,C,D):
    y = A * np.log(B*x + C) + D
    return y

real_OD = np.zeros((16,16))
real_OD[0,:] = [0.08,0.19,0.27,0.36,0.43,0.50,0.56,0.64,0.72,0.83,0.94,1.08,1.18,1.30,1.32,1.52]

for row in range(1,16):
    for col in range(0,16):
        real_OD[row,col] = real_OD[row-1,(col-1)%16]

detector_value = np.loadtxt(open("od_135_cal.txt", "rb"), delimiter=",", skiprows=0)

print(real_OD)
print(detector_value)
param0 = np.zeros([16,1])
param1 = np.zeros([16,1])
param2 = np.zeros([16,1])
param3 = np.zeros([16,1])

fig, axs = plt.subplots(4, 4)
for x in range(0,16):
    real_OD[:,x], detector_value[:,x] = zip(*sorted(zip(real_OD[:,x], detector_value[:,x])))

    axs[3-int(x/4), x%4].plot(real_OD[:,x],detector_value[:,x],'o')
    axs[3-int(x/4), x%4].set_title(("vial %d"%(x)))

    param, covar = curve_fit(fit_eqn9,real_OD[:,x], detector_value[:,x])

    x_test = np.linspace(0, 3, 300)

    fit_y = fit_eqn9(real_OD[:,x], param[0], param[1],param[2],param[3])

    axs[3-int(x/4), x%4].plot(real_OD[:,x], fit_y, '-')

    print(param)
    param0[x] = param[0]
    param1[x] = param[1]
    param2[x] = param[2]
    param3[x] = param[3]

filename = input("Please enter a name for this fit:\n")
f = open("od_cal.json",'w')
f.write("{\"name\": \"%s\", \"coefficients\": ["%filename)

for x in range (0,16):
    f.write("[%f,%f,%f,%f]"%(param0[x],param1[x],param2[x],param3[x]))
    if x < 15:
        f.write(",")

current_time = time.time()
        
f.write('], "type": "exp", "timeFit": %f, "active": true, "params": ["od_135"]}'%current_time)

f.close()

plt.show()
