"""
Will need to update the time parameter to reflect current time
"""

import numpy as np
import math
import matplotlib
import csv
import matplotlib.pyplot as plt
from random import *
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

def fit_eqn(x,A,B,C):
    y = A + B * np.log((C)/(x))
    return y

def fit_eqn2(x,A,B,C,D):
    y = A + B * np.log((C-D)/(x-D))
    return y

def fit_eqn3(x,A,B):
    y = A*x + B
    return y

def fit_eqn4(x,A,B,C):
    y = A*x*x + B*x + C
    return y

def fit_eqn5(x,A,B,C):
    y = A*10**(B*x) + C
    return y

def fit_eqn6(x,A,B,C,D):
    y = A*x**3 + B*x**2 + C*x + D
    return y

def fit_eqn7(x,A,B,C,D,E):
    y = A*x**4 + B*x**3 + C*x**2 + D*x + E
    return y

def fit_eqn32(c,x,y):
    return c[0] * x + c[1] - y

def fit_eqn8(x,A,B,C,D):
    y = A*10**(B*x + C) + D
    return y

def fit_eqn9(x,A,B,C,D):
    y = A * np.log(B*x + C) + D
    return y

real_OD = np.zeros((16,16))
real_OD[0,:] = [0.08,0.19,0.27,0.36,0.43,0.50,0.56,0.64,0.72,0.83,0.94,1.08,1.18,1.30,1.32,1.52]

for row in range(1,16):
    for col in range(0,16):
        real_OD[row,col] = real_OD[row-1,(col-1)%16]

detector_value = np.loadtxt(open("od_135_cal.txt", "rb"), delimiter=",", skiprows=0)
#detector_value = np.zeros((16,16))
#for row in range(0,16):
#    for col in range(0,16):
##            detector_value[row,col] = 17 + 2*(random() - 0.5)
#        else:
#            detector_value[row,col] = real_OD[row,col] + 2*(random() - 0.5)

#detector_value = real_OD

#real_OD2 = np.fliplr(real_OD)
#detector_value2 = np.flipud(detector_value)

#np.where(real_OD2==1,17,real_OD2)

print(real_OD)
print(detector_value)
param0 = np.zeros([16,1])
param1 = np.zeros([16,1])
param2 = np.zeros([16,1])
param3 = np.zeros([16,1])

fig, axs = plt.subplots(4, 4)
for x in range(0,16):
    real_OD[:,x], detector_value[:,x] = zip(*sorted(zip(real_OD[:,x], detector_value[:,x])))
    #axs[3-int(x/4), x%4].plot(detector_value[:,x],real_OD[:,x],'o')
    axs[3-int(x/4), x%4].plot(real_OD[:,x],detector_value[:,x],'o')
    axs[3-int(x/4), x%4].set_title(("vial %d"%(x)))
    #param, covar = curve_fit(fit_eqn9, detector_value[:,x],real_OD[:,x])
    param, covar = curve_fit(fit_eqn9,real_OD[:,x], detector_value[:,x])
    #param_robust = least_squares(fit_eqn32, [-0.0001, 3], loss='soft_l1', f_scale=0.1, args=(detector_value[:,x],real_OD[:,x]))
    x_test = np.linspace(0, 3, 300)

    #print(param_robust.x[0])
    #print(param_robust.x[1])

    #print(param_robust[0])y_robust = generate_data(x_test, *param_robust.x)

    #fit_y = fit_eqn8(detector_value[:,x], param[0], param[1],param[2],param[3])
    fit_y = fit_eqn9(real_OD[:,x], param[0], param[1],param[2],param[3])
    #y_robust = fit_eqn3(detector_value[:,x],param_robust.x[0],param_robust.x[1])

    axs[3-int(x/4), x%4].plot(real_OD[:,x], fit_y, '-')
    #axs[3-int(x/4), x%4].plot(detector_value[:,x], y_robust, '-')

    #param0[x] = param_robust.x[0]
    #param1[x] = param_robust.x[1]
    print(param)
    param0[x] = param[0]
    param1[x] = param[1]
    param2[x] = param[2]
    param3[x] = param[3]
    #axs[3-int(x/4), x%4].set_xlabel("y = %.2fx + %.2f" % (param[0], param[1]))

#plt.plot(real_OD2[:,0], detector_value2[:,0],'o')
#param, covar = curve_fit(fit_eqn4, real_OD2[:,0], detector_value2[:,0])
#print("y = %.2f + %.2f * log(%.2f-%.2f/x-%.2f)\n" %(param[0],param[1], param[2],param[3],param[3]))
#fit_y = fit_eqn4(real_OD2[:,0], param[0], param[1], param[2])
#plt.plot(real_OD2[:,0], fit_y, '-')
input = ("Please enter a name for this fit:\n")
f = open("od_cal.json",'w')
f.write("{\"name\": \"%s\", \"coefficients\": ["%input)

for x in range (0,16):
    f.write("[%f,%f,%f,%f]"%(param0[x],param1[x],param2[x],param3[x]))
    if x < 15:
        f.write(",")


f.write('], "type": "exp", "timeFit": 1634768074.8387222, "active": true, "params": ["od_135"]}')

f.close()

plt.show()
