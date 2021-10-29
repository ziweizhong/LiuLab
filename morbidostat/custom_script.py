#!/usr/bin/env python3
"""
2021-10-24 change log
- create a confirmation for pumps were run, and resends the message if they were not

2021-10-23 change log
- added initial power levels
- added additional operational mode for morbidostat
- added GLOBAL_VIALS parameter to set active vials

Possible things to update in the future:
- robust curve fitting (not working for non-linear curves)
"""
import numpy as np
import logging
import os.path
import time
import math
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime

# logger setup
logger = logging.getLogger(__name__)

##### USER DEFINED GENERAL SETTINGS #####

#set new name for each experiment, otherwise files will be overwritten
EXP_NAME = 'expt_ZZ_20211023_test_05'
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
OPERATION_MODE = 'morbidostat' # TURBIDOSTAT AND CHEMOSTAT FUNCTION HAS BEEN REMOVED FROM THIS VERSION
# if using a different mode, name your function as the OPERATION_MODE variable

GLOBAL_VIALS = [0,2]

##### END OF USER DEFINED GENERAL SETTINGS #####




def morbidostat(eVOLVER, input_data, vials, elapsed_time):
    # Timer
    #print(datetime.now().time())

    # Global paramemters for PID
    VOLUME =  30.0                              # Volume assumed for dilution events (mL)
    DILUTION_OD = 1                         # OD at which dilutions start occurring
    MAX_OD = 2.0                                # High OD safeguard to prevent culture from leaving linear range
    SUPER_OD = 3.0
    DOUBLING_TIME = 6.0                         # User defined target doubling time
    TARGET_GR = math.log(2.0)/DOUBLING_TIME     # Calculated based off the user defined doubling time
    STEPS_PER_HOUR = 4                          # User defined number of dilution events per hour. Volumes are adjusted accordingly
    DILUTION_RATE = float(VOLUME / DOUBLING_TIME) # mL / hour
    VOLUME_PER_DILUTION = float(DILUTION_RATE/ STEPS_PER_HOUR)
    TIME_BETWEEN_DILUTIONS = float(60/STEPS_PER_HOUR) #minutes
    TIME_OUT = VOLUME_PER_DILUTION + 3          # seconds to run pump for output

    # PID Parameters
    # Kp, Ki, Kd may need to be changed for different selections
    Kp =  [0.7] * len(vials)
    Ki =  [0.05] * len(vials)
    Kd =  [0.2] * len(vials)
    pid_offset = [0.0] * len(vials)


    current_OD_data = input_data['transformed']['od']
    save_path = os.path.dirname(os.path.realpath(__file__))

    mstat_vials = GLOBAL_VIALS #vials is all 16, can set to different range (ex. [0,1,2,3]) to only trigger tstat on those vials


    if elapsed_time > 0: # initial time for growth, hours

        #Load pump log (all recorded in vial 0 log across all 16, since pump time occurs simultaneously for all)
        file_path =  "%s/%s/pump_log/vial00_pump_log.txt" % (save_path,EXP_NAME)
        data = np.genfromtxt(file_path, delimiter=',')
        last_pump = data[len(data)-1][0]

        confirmed_pump_path =  "%s/%s/pump_log/confirmed_pump_log.txt" % (save_path,EXP_NAME)
        data = np.genfromtxt(confirmed_pump_path, delimiter=',')
        last_confirmed_pump = data[len(data)-1][0]

        # checks if dilution actually occured by ensuring OD measurements went down
        # If not, will resend dilution message
        if (last_pump - last_confirmed_pump)*3600 > 0 and (elapsed_time - last_pump)*3600 > 120:
            print("Last pump:%f"%last_pump)
            print("Last confirmed pump:%f"%last_confirmed_pump)
            message_file_path = "%s/%s/lastMessage.txt" % (save_path,EXP_NAME)
            message_file = open(message_file_path,'r')
            MESSAGE = message_file.read().strip().split(',')
            message_file.close()

            if MESSAGE == ['--'] * 48:
                #write current time as last confirmed pump
                confirmed_pump_file = open(confirmed_pump_path,'a+')
                confirmed_pump_file.write("%f,%f,%f\n"%(elapsed_time,od_before, od_after))
                confirmed_pump_file.close()
            else:
                # get vials that were diluted
                diluted_vials_path = "%s/%s/diluted_vials.txt" % (save_path,EXP_NAME)
                diluted_vials_file = open(diluted_vials_path,'r')
                diluted_vials = diluted_vials_file.read().strip().split(',')
                diluted_vials = list(filter(None,diluted_vials))
                diluted_vials_file.close()

                # obtains the first vial that underwent dilution
                od_path = "%s/%s/OD/vial%d_OD.txt" % (save_path,EXP_NAME,int(diluted_vials[0]))
                od_data = np.genfromtxt(od_path, delimiter=',')


                # finds the index of the last dilution
                time_data = od_data[:,0]
                idx = np.searchsorted(time_data,last_pump,side="left")

                # determine OD before the last dilutions
                od_before = np.mean(od_data[(idx-4):idx,1])
                # determine OD after the last dilution
                od_after = np.mean(od_data[(idx+2):(idx+5),1])
                print("OD before: %f" %od_before)
                print("OD after: %f" %od_after)
                # checks that OD went down after dilution
                if od_after / od_before < 0.98:
                    # dilution occurred
                    print("Dilution confirmed!")
                    # write current time as last confirmed pump
                    confirmed_pump_file = open(confirmed_pump_path,'a+')
                    confirmed_pump_file.write("%f,%f,%f\n"%(elapsed_time,od_before, od_after))
                    confirmed_pump_file.close()
                else:
                    # dilution did not occur
                    # resend message for dilutionsif MESSAGE != ['--'] * 48:
                    print("Dilution command resent!")
                    eVOLVER.fluid_command(MESSAGE)

                    # Updates pump file
                    pump_path =  "%s/%s/pump_log/vial00_pump_log.txt" % (save_path,EXP_NAME)
                    pump_file = open(pump_path,"a+")
                    pump_file.write("%f,%f\n" %  (elapsed_time,elapsed_time))
                    pump_file.close()

        # if enough time passed for the next dilution to occur
        if ((elapsed_time - last_pump)*3600) > (TIME_BETWEEN_DILUTIONS * 60): # convert to seconds and compare
            # determine which pumps to fire (e.g. drug or no drug) based on current OD
            #cleanCommand = 0;

            MESSAGE = ['--'] * 48
            diluted_vials = []

            for x in mstat_vials:
                #pumpCommand = 0
                #cleanCommand = cleanCommand + control[x+16]

                #pid_pump_command_no_drug = 0
                #pid_pump_command_drug   = 0

                # 1. Reads GR from last dilution time
                gr_path =  "%s/%s/growth_rate/vial%d_growth_rate.txt" % (save_path,EXP_NAME,x)
                gr_data = np.genfromtxt(gr_path, delimiter=',')
                last_gr = gr_data[len(gr_data)-1][1]

                # 2. Reads drug concentration from last dilution time
                drug_log_path =  "%s/%s/drug_log/vial%d_drug_log.txt" % (save_path,EXP_NAME,x)
                drug_log_data = np.genfromtxt(drug_log_path, delimiter=',')
                drug_conc = drug_log_data[len(drug_log_data)-1][1]

                # 3. Read and calculates the current OD (average over 5 measurements)
                od_path =  "%s/%s/OD/vial%d_OD.txt" % (save_path,EXP_NAME,x)
                od_data = np.genfromtxt(od_path, delimiter=',')

                # 4. Reads PID parameters from last round
                pid_path = "%s/%s/PIDLog/vial%d_PIDLog.txt" % (save_path,EXP_NAME,x)
                pid_data = np.genfromtxt(pid_path, delimiter=',')
                last_i_error = pid_data[len(pid_data)-1][2]
                last_p_error = pid_data[len(pid_data)-1][1]
                last_vial_pump = pid_data[len(pid_data)-1][0]

                # 5. Reads the offset parameter from last round
                offset_log_path = "%s/%s/offset_log/vial%d_offset.txt" % (save_path,EXP_NAME,x)
                offset_log_data = np.genfromtxt(offset_log_path, delimiter=',')
                if len(offset_log_data) > 2:
                    pid_offset[x] = offset_log_data[len(offset_log_data)-1][1]

                # CALCULATIONS #
                #Find index of last dilution event
                time_data = od_data[:,0]
                idx = np.searchsorted(time_data,last_pump,side="left")

                # chooses points after last dilution with a 5 tick offset to allow cultures to re-equilibrate
                od_window =   od_data[idx+5:,1]
                time_window = od_data[idx+5:,0]

                # resets time window to 0 for the first data points
                # makes the "A" constant in exp_func more meaningful
                time_window = time_window - time_window[0]

                # Removes NaNs and infs
                time_window = time_window[~np.isnan(od_window)]
                od_window = od_window[~np.isnan(od_window)]
                time_window = time_window[~np.isinf(od_window)]
                od_window = od_window[~np.isinf(od_window)]

                #print(time_window)
                #print(od_window)

                # Calculates OD for a vial (across 5 measurements)
                avg_od = 0
                for n in range(1,6):
                    avg_od = avg_od + (od_window[len(od_window)-n]/5)

                # Calculates OD after previous dilution (across 5 measurements)
                od_after_dilution = 0
                for n in range(0,3):
                    od_after_dilution = od_after_dilution + (od_window[n]/3)

                # Cleans up OD readings for growth rate calculation
                # This only occurs of OD < 0
                if od_after_dilution < 0 or avg_od < 0:
                    growth_rate = 0
                    gr_old = 0
                    od_fit_guess = 0
                else:
                    # WILL NEED TO CHANGE PUMP_LOG
                    """
                    file_name =  "vial{0}_curvefit_log.txt".format(x)
                    file_path = os.path.join(save_path, EXP_NAME, 'pump_log', file_name)
                    text_file = open(file_path, "a+")
                    text_file.write("\n\n" + str(elapsed_time) + "\n")
                    text_file.write(str(time_window))
                    text_file.write(str(od_window))
                    """
                    # Fits the OD and time to a exponential curve, func
                    #print(time_window)
                    #print(od_window)
                    param, converg = curve_fit(exp_func,time_window,od_window,bounds=([od_after_dilution - 0.05, -0.5],[od_after_dilution + 0.05, 0.5]))
                    param_robust = least_squares(exp_func_robust, [1, 0.1], loss='soft_l1', f_scale=0.1, args=(time_window,od_window))

                    #if x == 1:
                    #    plt.scatter(time_window, od_window)
                    #    plt.show()
                        #zipped = zip(time_window,od_window)
                        #np.savetxt('x.csv',time_window)
                        #np.savetxt('y.csv',od_window)

                    print("Original constant calculation: %.4f" %(param[0]))
                    print("Original growth rate calculation: %.4f" %(param[1]))
                    print("Robust constant calculation: %.4f" %(param_robust.x[0]))
                    print("Robust growth rate calculation: %.4f" %(param_robust.x[1]))

                    growth_rate = param[1]
                    #growth_rate = param.robust.x[1]

                # If the vial has been below the dilution OD for 4 rounds, and the culture has been growning for a while:
                # Dilute with no drug media to wash out drug.
                # Otherwise, dilute normally
                """
                #Removing this to save time/space #ZZ#
                if (elapsed_time - last_vial_pump) > (TIME_BETWEEN_DILUTIONS*4/60) and drug_conc > 0.05 and elapsed_time > 144:
                    print ("!!!!RESCUE DILUTION SHOULD HAVE OCCURRED!!!!")

                    # OLD PUMP COMMAND - #ZZ#: Ensure that new pump command is correct
                    #pumpCommand = pumpCommand + control[x] + control[x+16]
                    #MESSAGE   = "%s,0,%f," % ("{0:b}".format(pumpCommand) , 3)
                    #eVOLVER_module.fluid_command(MESSAGE, x, elapsedTime, TIME_BETWEEN_DILUTIONS *60, EXP_NAME, 3, 'n')

                    MESSAGE[x] = str(VOLUME_PER_DILUTION)
                    MESSAGE[x+16] == str(TIME_OUT)
                    # Logging pump commands
                    file_name =  "vial{0}_pump_log.txt".format(x)
                    file_path = os.path.join(save_path, EXP_NAME, 'pump_log', file_name)

                    text_file = open(file_path, "a+")
                    text_file.write("{0},{1},{2}\n".format(elapsed_time, VOLUME_PER_DILUTION,'0'))
                    text_file.close()

                    # Logging 0's for PID values, indicating a rescue dilution
                    pid_file = open(pid_path, "a+")
                    pid_file.write("%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (elapsed_time, 0, 0, 0, 0, 0, 0, 0, 0))
                    pid_file.close()

                    # Loggin drug concentration
                    drug_conc = (VOLUME * drug_conc)/(VOLUME + VOLUME_PER_DILUTION)
                    drug_file = open(drug_log_path,"a+")
                    drug_file.write("%f,%f\n" %  (elapsed_time, drug_conc))
                    drug_file.close()

                elif avg_od > DILUTION_OD:
                """
                if avg_od > DILUTION_OD:
                    # PID Calculations
                    p_err = 0
                    i_err = 0
                    d_err = 0

                    # pid_control is the % of dilution to be conducted with drug media
                    # will be calculated from errors and constants
                    pid_control = 0

                    # Calculates proportional error, e(t)
                    p_err = growth_rate - TARGET_GR

                    # Calculates Integrate[e(t), dt, t0, t] with trapezoid rule
                    # The integral is calculated every time the proportional
                    # error switches signs (so that every time that e[t] = 0, the integral term resets)
                    if p_err == 0 or last_p_error / p_err < 0:
                        i_err = (growth_rate + last_gr  - 2*TARGET_GR) / STEPS_PER_HOUR
                    else:
                        i_err = last_i_error + (growth_rate + last_gr  - 2*TARGET_GR) / STEPS_PER_HOUR

                    # Calculates d e(t)/dt
                    d_err = (growth_rate - last_gr) * STEPS_PER_HOUR

                    # If the machine has been acting normally, if the growth rate dips below the target growth rate,
                    # a new offset is calculated and used for future rounds of control
                    if elapsed_time > 72 and p_err < 0 and last_p_error > 0:
                        #print ("OFFSET LOG DIAGNOSTIC elapsed_time: %f" %(elapsed_time))
                        #print ("OFFSET LOG DIAGNOSTIC p_err: %f" %(p_err))
                        #print ("OFFSET LOG DIAGNOSTIC last_p_error: %f" %(last_p_error))
                        #print ("OFFSET LOG DIAGNOSTIC len(drug_log_data): %f" %(len(drug_log_data)))
                        new_pid_offset = 0
                        for n in range(6,11):
                            new_pid_offset = new_pid_offset + (drug_log_data[len(drug_log_data)-n][1])/5
                        pid_offset[x] = new_pid_offset
                        offset_file = open(offset_log_path, "a+")
                        offset_file.write("%f, %f\n" % (elapsed_time, new_pid_offset))
                        offset_file.close()

                    # Calculates the drug dilution percentage
                    pid_control = pid_offset[x] + p_err * Kp[x] + i_err * Ki[x] + d_err * Kd[x]

                    # Cleans up the controls in case they were calculated to be greater than 1 or less than 0
                    if pid_control > 1:
                        pid_control = 1
                    elif pid_control < 0:
                        pid_control = 0

                    # Option to dilute twice if the OD is > MAX_OD. This is to keep the culture within the linear range of the detectors
                    num_dils = 1
                    if avg_od > SUPER_OD:
                        num_dils = 8
                        #print ("Quad Dilution Occurred; num_dils = 4")
                    elif avg_od > MAX_OD:
                        num_dils = 2
                        #print ("Double Dilution Occurred; num_dils = 2")

                    # NEW PUMP Command
                    # No drug media influx
                    MESSAGE[x] = str(num_dils * VOLUME_PER_DILUTION * (1-pid_control))
                    # Drug media influx
                    MESSAGE[x+32] = str(num_dils * VOLUME_PER_DILUTION * (pid_control))
                    # Efflux pump
                    MESSAGE[x+16] = str(num_dils * VOLUME_PER_DILUTION + 5)

                    diluted_vials.append(str(x))

                    # Calculates the new drug concentration of the vial
                    #drug_conc = (VOLUME * drug_conc + pid_control * num_dils * VOLUME_PER_DILUTION/flow_rate)/(VOLUME + num_dils * VOLUME_PER_DILUTION/flow_rate)
                    drug_conc = (VOLUME * drug_conc + pid_control * num_dils * VOLUME_PER_DILUTION)/(VOLUME + num_dils * VOLUME_PER_DILUTION)

                    # Updates the PID parameter file
                    pid_file = open(pid_path, "a+")
                    pid_file.write("%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (elapsed_time, p_err, i_err, d_err, Kp[x] * p_err, Ki[x] * i_err, Kd[x] * d_err,pid_offset[x], pid_control))
                    pid_file.close()

                elif avg_od < DILUTION_OD:
                    p_err = 0
                    i_err = 0
                    d_err = 0
                    pid_offset[x] = 0
                    pid_control = 0
                    new_d_err = 0
                    gr_smooth = 0

                # 1. Update GR file with current GR
                GRFile = open(gr_path,"a+")
                GRFile.write("%f,%s\n" %  (elapsed_time, growth_rate))
                GRFile.close()

                # 2. Updates drug concentration file only if a dilution event happened
                if avg_od > DILUTION_OD:
                    drug_file = open(drug_log_path,"a+")
                    drug_file.write("%f,%f\n" %  (elapsed_time, drug_conc))
                    drug_file.close()

                # 3. Updates a log file with an assortment of values for plotting
                log_path = "%s/%s/logs/vial%d_log.txt" % (save_path,EXP_NAME,x)
                log_file = open(log_path, "a+")
                log_file.write ("Vial %d\n" %(x))
                log_file.write ("Average OD: %f\n" %(avg_od))
                log_file.write ("Target Growth Rate: %f\n" %(TARGET_GR))
                log_file.write ("Growth Rate: %f\n" %(growth_rate))
                log_file.write ("p_err: %f\n" %(p_err))
                log_file.write ("i error: %f\n" %(i_err))
                log_file.write ("d_err: %f\n" %(d_err))
                log_file.write ("PID Offset: %f\n" %(pid_offset[x]))
                log_file.write ("pid_control: %f\n" %(pid_control))
                log_file.close()


                # Terminal display values for each vial:
                print ("Vial %d" %(x))
                print ("Average OD: %f" %(avg_od))
                print ("Target Growth Rate: %f" %(TARGET_GR))
                print ("Growth Rate: %f" %(growth_rate))
                print ("last p error: %f" %(last_p_error))
                print ("p_err: %f" %(p_err))
                print ("i error: %f" %(i_err))
                print ("d_err: %f" %(d_err))
                print ("PID Offset: %f" %(pid_offset[x]))
                print ("pid_control: %f" %(pid_control))


                # Vial_log_file accounts for all Terminal output
                """
                log_path = "%s/%s/logs/vial00_log.txt" % (save_path,EXP_NAME)
                log_file = open(log_path, "a+")
                log_file.write ("\nElapsed Time: %f\n" %(elapsed_time))
                log_file.write ("Vial %d\n" %(x))
                log_file.write ("Average OD: %f\n" %(avg_od))
                log_file.write ("Target Growth Rate: %f\n" %(TARGET_GR))
                log_file.write ("Growth Rate: %f\n" %(growth_rate))
                log_file.write ("p_err: %f\n" %(p_err))
                log_file.write ("i error: %f\n" %(i_err))
                log_file.write ("d_err: %f\n" %(d_err))
                log_file.write ("PID Offset: %f\n" %(pid_offset[x]))
                log_file.write ("pid_control: %f\n" %(pid_control))
                #log_file.write ("Message Drug: %s\n" %(MESSAGEDrug))
                #log_file.write ("Message NoDrug: %s\n\n" %(MESSAGENoDrug))
                log_file.close()
                """

            print(MESSAGE[0:4])
            print(MESSAGE[16:20])
            print(MESSAGE[32:36])

            message_path =  "%s/%s/lastMessage.txt" % (save_path,EXP_NAME)
            text_file = open(message_path,"w")
            for index in range(0,48):
                if index != 47:
                    text_file.write(MESSAGE[index]+',')
                else:
                    text_file.write(MESSAGE[index])
            text_file.close()

            diluted_vials_path =  "%s/%s/diluted_vials.txt" % (save_path,EXP_NAME)
            text_file = open(diluted_vials_path,"w")
            for item in diluted_vials:
                text_file.write(item+',')
            text_file.close()

            if MESSAGE != ['--'] * 48:
                print("Dilution performed!")
                eVOLVER.fluid_command(MESSAGE)


                # Updates pump file
            pump_path =  "%s/%s/pump_log/vial00_pump_log.txt" % (save_path,EXP_NAME)
            pump_file = open(pump_path,"a+")
            pump_file.write("%f,%f\n" %  (elapsed_time,elapsed_time))
            pump_file.close()



    #print(datetime.now().time())

def turbidostat(eVOLVER, input_data, vials, elapsed_time):
    OD_data = input_data['transformed']['od']

    ##### USER DEFINED VARIABLES #####

    turbidostat_vials = GLOBAL_VIALS #vials is all 16, can set to different range (ex. [0,1,2,3]) to only trigger tstat on those vials
    stop_after_n_curves = np.inf #set to np.inf to never stop, or integer value to stop diluting after certain number of growth curves
    OD_values_to_average = 6  # Number of values to calculate the OD average

    lower_thresh = [1] * len(vials) #to set all vials to the same value, creates 16-value list
    upper_thresh = [4] * len(vials) #to set all vials to the same value, creates 16-value list

    #Alternatively, use 16 value list to set different thresholds, use 9999 for vials not being used
    #lower_thresh = [0.2, 0.2, 0.3, 0.3, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999]
    #upper_thresh = [0.4, 0.4, 0.4, 0.4, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999]


    ##### END OF USER DEFINED VARIABLES #####


    ##### Turbidostat Settings #####
    #Tunable settings for overflow protection, pump scheduling etc. Unlikely to change between expts

    time_out = 5 #(sec) additional amount of time to run efflux pump
    pump_wait = 3 # (min) minimum amount of time to wait between pump events

    ##### End of Turbidostat Settings #####

    save_path = os.path.dirname(os.path.realpath(__file__)) #save path
    flow_rate = eVOLVER.get_flow_rate() #read from calibration file


    ##### Turbidostat Control Code Below #####

    # fluidic message: initialized so that no change is sent
    MESSAGE = ['--'] * 48
    for x in turbidostat_vials: #main loop through each vial

        # Update turbidostat configuration files for each vial
        # initialize OD and find OD path

        file_name =  "vial{0}_ODset.txt".format(x)
        ODset_path = os.path.join(save_path, EXP_NAME, 'ODset', file_name)
        data = np.genfromtxt(ODset_path, delimiter=',')
        ODset = data[len(data)-1][1]
        ODsettime = data[len(data)-1][0]
        num_curves=len(data)/2;

        file_name =  "vial{0}_OD.txt".format(x)
        OD_path = os.path.join(save_path, EXP_NAME, 'OD', file_name)
        data = eVOLVER.tail_to_np(OD_path, OD_values_to_average)
        average_OD = 0

        # Determine whether turbidostat dilutions are needed
        #enough_od_data = (len(data) > 7) #logical, checks to see if enough data points (couple minutes) for sliding window
        collecting_more_curves = (num_curves <= (stop_after_n_curves + 2)) #logical, checks to see if enough growth curves have happened

        if data.size != 0:
            # Take median to avoid outlier
            od_values_from_file = data[:,1]
            average_OD = float(np.median(od_values_from_file))

            #if recently exceeded upper threshold, note end of growth curve in ODset, allow dilutions to occur and growth_rate to be measured
            if (average_OD > upper_thresh[x]) and (ODset != lower_thresh[x]):
                text_file = open(ODset_path, "a+")
                text_file.write("{0},{1}\n".format(elapsed_time,
                                                   lower_thresh[x]))
                text_file.close()
                ODset = lower_thresh[x]
                # calculate growth rate
                eVOLVER.calc_growth_rate(x, ODsettime, elapsed_time)

            #if have approx. reached lower threshold, note start of growth curve in ODset
            if (average_OD < (lower_thresh[x] + (upper_thresh[x] - lower_thresh[x]) / 3)) and (ODset != upper_thresh[x]):
                text_file = open(ODset_path, "a+")
                text_file.write("{0},{1}\n".format(elapsed_time, upper_thresh[x]))
                text_file.close()
                ODset = upper_thresh[x]

            #if need to dilute to lower threshold, then calculate amount of time to pump
            if average_OD > ODset and collecting_more_curves:

                time_in = - (np.log(lower_thresh[x]/average_OD)*VOLUME)/flow_rate[x]

                if time_in > 20:
                    time_in = 20

                time_in = round(time_in, 2)

                file_name =  "vial{0}_pump_log.txt".format(x)
                file_path = os.path.join(save_path, EXP_NAME,
                                         'pump_log', file_name)
                data = np.genfromtxt(file_path, delimiter=',')
                last_pump = data[len(data)-1][0]
                if ((elapsed_time - last_pump)*60) >= pump_wait: # if sufficient time since last pump, send command to Arduino
                    logger.info('turbidostat dilution for vial %d' % x)
                    # influx pump
                    MESSAGE[x] = str(time_in)
                    # efflux pump
                    MESSAGE[x + 16] = str(time_in + time_out)

                    file_name =  "vial{0}_pump_log.txt".format(x)
                    file_path = os.path.join(save_path, EXP_NAME, 'pump_log', file_name)

                    text_file = open(file_path, "a+")
                    text_file.write("{0},{1}\n".format(elapsed_time, time_in))
                    text_file.close()
        else:
            logger.debug('not enough OD measurements for vial %d' % x)

    # send fluidic command only if we are actually turning on any of the pumps
    if MESSAGE != ['--'] * 48:
        print("Dilution occurred!!")
        print(MESSAGE[0:16])
        print(MESSAGE[16:32])
        print(MESSAGE[32:48])
        eVOLVER.fluid_command(MESSAGE)

        # your_FB_function_here() #good spot to call feedback functions for dynamic temperature, stirring, etc for ind. vials
    # your_function_here() #good spot to call non-feedback functions for dynamic temperature, stirring, etc.

    # end of turbidostat() fxn

def chemostat(eVOLVER, input_data, vials, elapsed_time):
    OD_data = input_data['transformed']['od_90']

    ##### USER DEFINED VARIABLES #####
    start_OD = 0 # ~OD600, set to 0 to start chemostate dilutions at any positive OD
    start_time = 0 #hours, set 0 to start immediately
    # Note that script uses AND logic, so both start time and start OD must be surpassed

    OD_values_to_average = 6  # Number of values to calculate the OD average
    chemostat_vials = GLOBAL_VIALS #vials is all 16, can set to different range (ex. [0,1,2,3]) to only trigger tstat on those vials

    rate_config = [0.5] * 16 #to set all vials to the same value, creates 16-value list
    #UNITS of 1/hr, NOT mL/hr, rate = flowrate/volume, so dilution rate ~ growth rate, set to 0 for unused vials

    #Alternatively, use 16 value list to set different rates, use 0 for vials not being used
    #rate_config = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6]

    ##### END OF USER DEFINED VARIABLES #####


    ##### Chemostat Settings #####

    #Tunable settings for bolus, etc. Unlikely to change between expts
    bolus = 0.5 #mL, can be changed with great caution, 0.2 is absolute minimum

    ##### End of Chemostat Settings #####

    save_path = os.path.dirname(os.path.realpath(__file__)) #save path
    flow_rate = eVOLVER.get_flow_rate() #read from calibration file
    period_config = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] #initialize array
    bolus_in_s = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] #initialize array


    ##### Chemostat Control Code Below #####

    for x in chemostat_vials: #main loop through each vial

        # Update chemostat configuration files for each vial

        #initialize OD and find OD path
        file_name =  "vial{0}_OD.txt".format(x)
        OD_path = os.path.join(save_path, EXP_NAME, 'OD', file_name)
        data = eVOLVER.tail_to_np(OD_path, OD_values_to_average)
        average_OD = 0
        #enough_od_path = (len(data) > 7) #logical, checks to see if enough data points (couple minutes) for sliding window

        if data.size != 0: #waits for seven OD measurements (couple minutes) for sliding window

            #calculate median OD
            od_values_from_file = data[:,1]
            average_OD = float(np.median(od_values_from_file))

            # set chemostat config path and pull current state from file
            file_name =  "vial{0}_chemo_config.txt".format(x)
            chemoconfig_path = os.path.join(save_path, EXP_NAME,
                                            'chemo_config', file_name)
            chemo_config = np.genfromtxt(chemoconfig_path, delimiter=',')
            last_chemoset = chemo_config[len(chemo_config)-1][0] #should t=0 initially, changes each time a new command is written to file
            last_chemophase = chemo_config[len(chemo_config)-1][1] #should be zero initially, changes each time a new command is written to file
            last_chemorate = chemo_config[len(chemo_config)-1][2] #should be 0 initially, then period in seconds after new commands are sent

            # once start time has passed and culture hits start OD, if no command has been written, write new chemostat command to file
            if ((elapsed_time > start_time) & (average_OD > start_OD)):

                #calculate time needed to pump bolus for each pump
                bolus_in_s[x] = bolus/flow_rate[x]

                # calculate the period (i.e. frequency of dilution events) based on user specified growth rate and bolus size
                if rate_config[x] > 0:
                    period_config[x] = (3600*bolus)/((rate_config[x])*VOLUME) #scale dilution rate by bolus size and volume
                else: # if no dilutions needed, then just loops with no dilutions
                    period_config[x] = 0

                if  (last_chemorate != period_config[x]):
                    print('Chemostat updated in vial {0}'.format(x))
                    logger.info('chemostat initiated for vial %d, period %.2f'
                                % (x, period_config[x]))
                    # writes command to chemo_config file, for storage
                    text_file = open(chemoconfig_path, "a+")
                    text_file.write("{0},{1},{2}\n".format(elapsed_time,
                                                           (last_chemophase+1),
                                                           period_config[x])) #note that this changes chemophase
                    text_file.close()
        else:
            logger.debug('not enough OD measurements for vial %d' % x)

        # your_FB_function_here() #good spot to call feedback functions for dynamic temperature, stirring, etc for ind. vials
    # your_function_here() #good spot to call non-feedback functions for dynamic temperature, stirring, etc.

    eVOLVER.update_chemo(input_data, chemostat_vials, bolus_in_s, period_config) #compares computed chemostat config to the remote one
    # end of chemostat() fxn


# def your_function_here(): # good spot to define modular functions for dynamics or feedback
def exp_func(x, a, b):
    # Exponential function for fit
    return a * np.exp(b * x)

def exp_func_robust(param,x,y):
    # exponential function for robust fit
    return param[0] * np.exp(param[1] * x) - y

if __name__ == '__main__':
    print('Please run eVOLVER.py instead')
