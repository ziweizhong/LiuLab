#!/usr/bin/env python3

import numpy as np
import logging
import os.path
import time
import math
from scipy.optimize import curve_fit
#from scipy.optimize import least_squares
#import matplotlib
#import matplotlib.pyplot as plt
from datetime import datetime
# logger setup
logger = logging.getLogger(__name__)

##### USER DEFINED GENERAL SETTINGS #####

#set new name for each experiment, otherwise files will be overwritten
EXP_NAME = 'expt_ZZ_20211021_OD_calibration'
EVOLVER_IP = '10.0.0.3'
EVOLVER_PORT = 8081

##### Identify pump calibration files, define initial values for temperature, stirring, volume, power settings

TEMP_INITIAL = [30] * 16 #degrees C, makes 16-value list
#Alternatively enter 16-value list to set different values
#TEMP_INITIAL = [30,30,30,30,32,32,32,32,34,34,34,34,36,36,36,36]

STIR_INITIAL = [6,9,6,6] * 4 #try 8,10,12 etc; makes 16-value list
#Alternatively enter 16-value list to set different values
#STIR_INITIAL = [7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10]

POWER_INITIAL = [4095] * 16

VOLUME =  25 #mL, determined by vial cap straw length
PUMP_CAL_FILE = 'pump_cal.txt' #tab delimited, mL/s with 16 influx pumps on first row, etc.
OPERATION_MODE = 'calibrate_OD' # TURBIDOSTAT AND CHEMOSTAT FUNCTION HAS BEEN REMOVED FROM THIS VERSION
# if using a different mode, name your function as the OPERATION_MODE variable

GLOBAL_VIALS = range(0,16)

##### END OF USER DEFINED GENERAL SETTINGS #####

def calibrate_OD(eVOLVER, input_data, vials, elapsed_time):
    # separates the data from the 90 deg and 135 deg detectors
    od_90_data = input_data['data'].get('od_90', None)
    od_135_data = input_data['data'].get('od_135', None)

    # gets the number of calibrations done so far by looking at the length of the file
    old_data = np.genfromtxt('od_90_cal.txt',delimiter=',')
    
    # if the number of calibrations saved is less than the total of 16, obtain more
    if old_data.ndim < 2 or old_data.shape[0] < 16:
        
        # outputs the raw data and gives the option to save
        print('OD 90 Data:')
        print(od_90_data)
        print('OD 135 Data:')
        print(od_135_data)
        choice = input("Data obtained! Do you want to save this data? (y/n)")
        
        if choice == 'y':
            
            # writes data to file
            with open('od_90_cal.txt','a+') as file:
                for item in od_90_data:
                    file.write(item + ',')
                file.write("\n")
            with open('od_135_cal.txt','a+') as file:
                for item in od_135_data:
                    file.write(item + ',')
                file.write("\n")
                
            # outputs the number of datapoints gathered so far
            if old_data.ndim < 2 or old_data.shape[0] < 16:
                if old_data.ndim == 0:
                    counter = 0
                elif old_data.ndim == 1:
                    counter = 1
                else:
                    counter = old_data.shape[0]
                print("\n\n___________________________________________________________")
                print ("Data %d of 16 saved!"% counter+1)
                if counter < 15:
                    print("Please advance all vails by 1 and move last vial to first vial")
        else:
            print("Data not saved! Please reobtain the data!")
            
        # If the last data point has not been gathered, give directions to advance vials
        if old_data.ndim < 2 or old_data.shape[0] < 15:
            input("Press enter when ready...")
            print("Waiting for data...")
    else:
        print("All calibrations obtained! Press Ctrl-C to quit!")


if __name__ == '__main__':
    print('Please run eVOLVER.py instead')
