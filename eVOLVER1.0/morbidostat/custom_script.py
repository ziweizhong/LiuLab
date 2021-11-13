import eVOLVER_module
import numpy as np
import os.path
import math
from scipy.optimize import curve_fit

###### Description of Code ##########
# This code is designed for eVOLVER to dynamically tune selection by altering the stringency of selection (AKA drug concentration)
# via diluting with a mixture of two choices of media; one media that is maximum drug concentration and one media that has no drug
# or a drug concentration at which selection first begins.
# Every 1 hour (or user definable timeline), a dilution event occurs to dilute a fixed amount based on the culture volume (an estimate),
# and the target growth rate (user defined).
"""
# PUMP TABLE
# PUMPS   VIAL    ROLE
# 0   0   NO DRUG
# 1   1   NO DRUG
# 2   2   NO DRUG
# 3   3   NO DRUG
# 4   4   NO DRUG
# 5   5   NO DRUG
# 6   6   NO DRUG
# 7   7   NO DRUG
# 8   0   DRUG
# 9   1   DRUG
# 10  2   DRUG
# 11  3   DRUG
# 12  4   DRUG
# 13  5   DRUG
# 14  6   DRUG
# 15  7   DRUG
# 16  0   WASTE
# 17  1   WASTE
# 18  2   WASTE
# 19  3   WASTE
# 20  4   WASTE
# 21  5   WASTE
# 22  6   WASTE
# 23  7   WASTE
# 24  12    NO DRUG
# 25  13    NO DRUG
# 26  12    DRUG
# 27  13    DRUG
# 28  12    WASTE
# 29  13    WASTE
# 30  N/A N/A
# 31  N/A N/A
"""
GLOBAL_VIALS = [0,1,2,3]

# def morbidostat (ODData, tempData, vials, elapsed_time, exp_name):
def morbidostat (current_OD_data, temp_data, vials, elapsed_time, exp_name):
    # Set Stir rates
    MESSAGE = "10,10,10,10   ,10,10,10,10,  10,10,10,10,  10,10,10,10,"
    eVOLVER_module.stir_rate(MESSAGE)

    control = np.power(2,range(0,32))
    # Calibrated flow rates
    flow_rate = 1.0 #ml/sec

    # Global paramemters for PID
    VOLUME =  30.0                              # Volume assumed for dilution events (mL)
    DILUTION_OD = 0.80                          # OD at which dilutions start occurring
    MAX_OD = 1.4                               # High OD safeguard to prevent culture from leaving linear range
    DOUBLING_TIME = 8.0                         # User defined target doubling time
    TARGET_GR = math.log(2.0)/DOUBLING_TIME     # Calculated based off the user defined doubling time
    STEPS_PER_HOUR = 2                          # User defined number of dilution events per hour. Volumes are adjusted accordingly
    DILUTION_RATE = float(VOLUME / DOUBLING_TIME) # mL / hour
    VOLUME_PER_DILUTION = float(DILUTION_RATE/ STEPS_PER_HOUR)
    TIME_BETWEEN_DILUTIONS = float(60/STEPS_PER_HOUR) #minutes

    # PID Parameters
    # Kp, Ki, Kd may need to be changed for different selections
    #Kp =         [0.07, 0.07,        0.07, 0.07,         0.07, 0.07,          0.07, 0.07]
    #Ki =         [0.05, 0.05,        0.05, 0.05,     0.05, 0.05,      0.05, 0.05]
    #Kd =         [0.2, 0.2,        0.2, 0.2,         0.2, 0.2,          0.2, 0.2]
    #pid_offset = [0.0, 0.0,    0.0, 0.0,     0.0, 0.0,      0.0, 0.0]
    Kp = [0.7] * 16
    Ki = [0.05] * 16
    Kd = [0.2] * 16
    pid_offset = [0.0] * 16

    
    if elapsed_time > 24: # initial time for growth, hours

        #Load pump log (all recorded in vial 0 log across all 16, since pump time occurs simultaneously for all)
        save_path = os.path.dirname(os.path.realpath(__file__))
        file_path =  "%s/%s/pump_log/vial0_pump_log.txt" % (save_path,exp_name)
        data = np.genfromtxt(file_path, delimiter=',')
        last_pump = data[len(data)-1][0]

        #Load pump log (all recorded in vial 0 log across all 16, since pump time occurs simultaneously for all)
        file_path =  "%s/%s/pump_log/vial00_pump_log.txt" % (save_path,exp_name)
        data = np.genfromtxt(file_path, delimiter=',')
        last_pump = data[len(data)-1][0]

        confirmed_pump_path =  "%s/%s/pump_log/confirmed_pump_log.txt" % (save_path,exp_name)
        data = np.genfromtxt(confirmed_pump_path, delimiter=',')
        last_confirmed_pump = data[len(data)-1][0]

        od_path = "%s/%s/OD/vial%d_OD.txt" % (save_path,exp_name,int(GLOBAL_VIALS[0]))
        od_data = np.genfromtxt(od_path, delimiter=',')
        last_data_time = od_data[len(od_data)-1][0]

        # To enable dilution check, will have to determine how to resent individual pump commands
        """
        # checks if dilution actually occured by ensuring OD measurements went down
        # If not, will resend dilution message
        if (last_pump - last_confirmed_pump)*3600 > 0 and (elapsed_time - last_pump)*3600 > 120:
            print("Last pump:%f"%last_pump)
            print("Last confirmed pump:%f"%last_confirmed_pump)
            message_file_path = "%s/%s/lastMessage.txt" % (save_path,exp_name)
            message_file = open(message_file_path,'r')
            MESSAGE = message_file.read().strip().split(',')
            message_file.close()

            # get vials that were diluted
            diluted_vials_path = "%s/%s/diluted_vials.txt" % (save_path,exp_name)
            diluted_vials_file = open(diluted_vials_path,'r')
            diluted_vials = diluted_vials_file.read().strip().split(',')
            diluted_vials = list(filter(None,diluted_vials))
            diluted_vials_file.close()

            # obtains the first vial that underwent dilution
            if diluted_vials != []:
                od_path = "%s/%s/OD/vial%d_OD.txt" % (save_path,exp_name,int(diluted_vials[0]))
                od_data = np.genfromtxt(od_path, delimiter=',')
            else:
                od_path = "%s/%s/OD/vial%d_OD.txt" % (save_path,exp_name,int(GLOBAL_VIALS[0]))
                od_data = np.genfromtxt(od_path, delimiter=',')

            # finds the index of the last dilution
            time_data = od_data[:,0]
            idx = np.searchsorted(time_data,last_pump,side="left")

            # determine OD before the last dilutions
            od_before = np.mean(od_data[(idx-4):idx,1])
            # determine OD after the last dilution
            od_after = np.mean(od_data[(idx+2):(idx+5),1])

            if MESSAGE == ['--'] * 48:
                #write current time as last confirmed pump
                confirmed_pump_file = open(confirmed_pump_path,'a+')
                confirmed_pump_file.write("%f,%f,%f\n"%(elapsed_time,od_before, od_after))
                confirmed_pump_file.close()
            else:

                print("OD before: %f" %od_before)
                print("OD after: %f" %od_after)
                # checks that OD went down after dilution
                if od_after / od_before < 0.995:
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
                    pump_path =  "%s/%s/pump_log/vial00_pump_log.txt" % (save_path,exp_name)
                    pump_file = open(pump_path,"a+")
                    pump_file.write("%f,%f\n" %  (elapsed_time,elapsed_time))
                    pump_file.close()
        """

        # if enough time passed for the next dilution to occur
        if ((elapsed_time- last_pump)*3600) > (TIME_BETWEEN_DILUTIONS * 60) and (last_data_time - last_pump)*3600 > 300: # convert to seconds and compare

            # determine which pumps to fire (e.g. drug or no drug) based on current OD
            clean_command = 0

            for x in GLOBAL_VIALS:
                pump_command = 0
                clean_command = clean_command + control[x+16]

                pid_pump_command_no_drug = 0
                pid_pump_command_drug   = 0

                # 1. Reads GR from last dilution time
                gr_path =  "%s/%s/growth_rate/vial%d_growth_rate.txt" % (save_path,exp_name,x)
                gr_data = np.genfromtxt(gr_path, delimiter=',')
                last_gr = gr_data[len(gr_data)-1][1]

                # 2. Reads drug concentration from last dilution time
                drug_log_path =  "%s/%s/drug_log/vial%d_drug_log.txt" % (save_path,exp_name,x)
                drug_log_data = np.genfromtxt(drug_log_path, delimiter=',')
                drug_conc = drug_log_data[len(drug_log_data)-1][1]

                # 3. Read and calculates the current OD (average over 5 measurements)
                od_path =  "%s/%s/OD/vial%d_OD.txt" % (save_path,exp_name,x)
                od_data = np.genfromtxt(od_path, delimiter=',')

                # 4. Reads PID parameters from last round
                pid_path = "%s/%s/PIDLog/vial%d_PIDLog.txt" % (save_path,exp_name,x)
                pid_data = np.genfromtxt(pid_path, delimiter=',')
                last_i_error = pid_data[len(pid_data)-1][2]
                last_p_error = pid_data[len(pid_data)-1][1]
                last_vial_pump = pid_data[len(pid_data)-1][0]

                # 5. Reads the offset parameter from last round
                offset_log_path = "%s/%s/offset_log/vial%d_offset.txt" % (save_path,exp_name,x)
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

                    # Fits the OD and time to a exponential curve, func
                    param, converg = curve_fit(exp_func,time_window,od_window,bounds=([od_after_dilution - 0.05, 0],[od_after_dilution + 0.05, 0.5]))

                    growth_rate = param[1]

                # If the vial has been below the dilution OD for 4 rounds, and the culture has been growning for a while:
                # Dilute with no drug media to wash out drug.
                # Otherwise, dilute normally
                """
                Removing this for now. Can add in later.
                if (elapsed_time - lastVialPump) > (TIME_BETWEEN_DILUTIONS*4/60) and drugConc > 0.05 and elapsed_time > 144:
                    print ("!!!!RESCUE DILUTION SHOULD HAVE OCCURRED!!!!")
                    pump_command = pump_command + control[x] + control[x+16]
                    MESSAGE   = "%s,0,%f," % ("{0:b}".format(pump_command) , 3)
                    eVOLVER_module.fluid_command(MESSAGE, x, elapsed_time, TIME_BETWEEN_DILUTIONS *60, exp_name, 3, 'n')
                    pid_file = open(pid_path, "a+")
                    pid_file.write("%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (elapsed_time, 0, 0, 0, 0, 0, 0, 0, 0))
                    pid_file.close()

                    drugConc = (VOLUME * drugConc)/(VOLUME + 3)
                    drug_file = open(drug_log_path,"a+")
                    drug_file.write("%f,%f\n" %  (elapsed_time, drugConc))
                    drug_file.close()
                elif avg_od > DILUTION_OD:
                """

                # PID Calculations
                p_err = 0
                i_err = 0
                d_err = 0

                MESSAGE_drug = "0,0,0"
                MESSAGE_no_drug = "0,0,0"

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
                # Calculates d e(t)/dt, but with a smoothed fit
                # Smoothing using Savitzky-Golay filter, using 5 points
                # Also smoothes the growth rate
                # However, this means that the smoothed function is two time points behind
                # Can/should this still be used in the PID? #ZZ#
                if len(gr_data) > 4:
                    new_d_err = STEPS_PER_HOUR / 12 * (gr_data[len(gr_data)-5][1] - 8 * gr_data[len(gr_data)-4][1]  + 8 * gr_data[len(gr_data)-2][1] - gr_data[len(gr_data)-1][1])

                    gr_smooth = 1 / 35 * (12 * gr_data[len(gr_data)-4][1] + 17 * gr_data[len(gr_data)-3][1] + 12 * gr_data[len(gr_data)-2][1] - 3 * gr_data[len(gr_data)-1][1] - 3 * gr_data[len(gr_data)-5][1])

                    d_err_smooth_window = [gr_data[len(gr_data)-4][1],gr_data[len(gr_data)-3][1],gr_data[len(gr_data)-2][1],gr_data[len(gr_data)-1][1]]
                    d_err_time_window = np.arange(0.0,(4)/STEPS_PER_HOUR,float(1/STEPS_PER_HOUR))

                    derr_param, derr_converg = curve_fit(lin_func,d_err_time_window,d_err_smooth_window)
                    d_err_smooth = derr_param[0]

                else:
                    new_d_err = 0
                    gr_smooth = 0

                
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

                if avg_od > DILUTION_OD:
                    # PID Pump Command
                    pid_pump_command_no_drug = pid_pump_command_no_drug + control[x] + control[x+16]
                    pid_pump_command_drug   = pid_pump_command_drug   + control[x+8] + control[x+16]

                    # Option to dilute twice if the OD is > MAX_OD. This is to keep the culture within the linear range of the detectors
                    num_dils = 1
                    if avg_od > MAX_OD:
                        num_dils = 2
                        print ("Double Dilution Occurred; num_dils = 2")

                    # Sends dilution command to RaspberryPi
                    MESSAGE_no_drug = "%s,0,%f," % ("{0:b}".format(pid_pump_command_no_drug) , num_dils * VOLUME_PER_DILUTION/flow_rate * ( 1 - pid_control))
                    MESSAGE_drug   = "%s,0,%f," % ("{0:b}".format(pid_pump_command_drug) , num_dils * VOLUME_PER_DILUTION/flow_rate * pid_control)
                    if pid_control < 1:
                        eVOLVER_module.fluid_command(MESSAGE_no_drug, x, elapsed_time, TIME_BETWEEN_DILUTIONS *60, exp_name, num_dils * VOLUME_PER_DILUTION/flow_rate * ( 1 - pid_control), 'n')
                        print ("Message NoDrug: %s" %(MESSAGE_no_drug))
                    if pid_control > 0:
                        eVOLVER_module.fluid_command(MESSAGE_drug, x, elapsed_time, TIME_BETWEEN_DILUTIONS *60, exp_name, num_dils * VOLUME_PER_DILUTION/flow_rate * pid_control, 'n')
                        print ("Message Drug: %s" %(MESSAGE_drug))

                    # Calculates the new drug concentration of the vial
                    drugConc = (VOLUME * drugConc + pid_control * num_dils * VOLUME_PER_DILUTION/flow_rate)/(VOLUME + num_dils * VOLUME_PER_DILUTION/flow_rate)
                else:
                    pid_control = 0

                # Updates the PID parameter file
                pid_file = open(pid_path, "a+")
                pid_file.write("%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (elapsed_time, p_err, i_err, d_err, Kp[x] * p_err, Ki[x] * i_err, Kd[x] * d_err,pid_offset[x], pid_control))
                pid_file.close()


                # 1. Update GR file with current GR
                gr_file = open(gr_path,"a+")
                gr_file.write("%f,%s\n" %  (elapsed_time, growth_rate))
                gr_file.close()

                # 2. Updates drug concentration file only if a dilution event happened
                if avg_od > DILUTION_OD:
                    drug_file = open(drug_log_path,"a+")
                    drug_file.write("%f,%f\n" %  (elapsed_time, drugConc))
                    drug_file.close()

                # 3. Updates a log file with an assortment of values for plotting
                log_path = "%s/%s/logs/vial%d_log.txt" % (save_path,exp_name,x)
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
                log_file.write ("Message Drug: %s\n" %(MESSAGE_drug))
                log_file.write ("Message NoDrug: %s\n\n" %(MESSAGE_no_drug))
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

                # Vial13_log_file accounts for all Terminal output
                log_path = "%s/%s/logs/vial13_log.txt" % (save_path,exp_name)
                log_file = open(log_path, "a+")
                log_file.write ("Elapsed Time: %f" %(elapsed_time))
                log_file.write ("Vial %d\n" %(x))
                log_file.write ("Average OD: %f\n" %(avg_od))
                log_file.write ("Target Growth Rate: %f\n" %(TARGET_GR))
                log_file.write ("Growth Rate: %f\n" %(growth_rate))
                log_file.write ("p_err: %f\n" %(p_err))
                log_file.write ("i error: %f\n" %(i_err))
                log_file.write ("d_err: %f\n" %(d_err))
                log_file.write ("PID Offset: %f\n" %(pid_offset[x]))
                log_file.write ("pid_control: %f\n" %(pid_control))
                log_file.write ("Message Drug: %s\n" %(MESSAGE_drug))
                log_file.write ("Message NoDrug: %s\n\n" %(MESSAGE_no_drug))
                log_file.close()

            # Runs output pumps for 5 seonds to ensure volume of vial is constant
            CLEAN_MESSAGE = "%s,0,5," % ("{0:b}".format(clean_command))
            print ("CLEAN_MESSAGE: %s" %(CLEAN_MESSAGE))
            eVOLVER_module.fluid_command(CLEAN_MESSAGE, 0, elapsed_time, TIME_BETWEEN_DILUTIONS *60, exp_name, 5,'y')

            
            #Code to complement confirm last pump.
            """
            message_path =  "%s/%s/lastMessage.txt" % (save_path,exp_name)
            text_file = open(message_path,"w")
            for index in range(0,48):
                if index != 47:
                    text_file.write(MESSAGE[index]+',')
                else:
                    text_file.write(MESSAGE[index])
            text_file.close()

            diluted_vials_path =  "%s/%s/diluted_vials.txt" % (save_path,exp_name)
            text_file = open(diluted_vials_path,"w")
            for item in diluted_vials:
                text_file.write(item+',')
            text_file.close()

            # Updates pump file
            pump_path =  "%s/%s/pump_log/vial00_pump_log.txt" % (save_path,exp_name)
            pump_file = open(pump_path,"a+")
            pump_file.write("%f,%f\n" %  (elapsed_time,elapsed_time))
            pump_file.close()
            """


def exp_func(x, a, b):
    # Exponential function for fit
    return a * np.exp(b * x)

def lin_func(x,a,b):
    return a * x + b

def measure_growth_rate(x, save_path, exp_name):
    ## Load files from proper directories
    ODPath =  "%s/%s/OD/vial%d_OD.txt" % (save_path,exp_name,x)
    pump_path =  "%s/%s/pump_log/vial0_pump_log.txt" % (save_path,exp_name)
    ODData = np.genfromtxt(ODPath, delimiter=',')
    pump_data = np.genfromtxt(pump_path, delimiter=',')

    ### Gets time frame from last diution event to current event
    time_buffer = 3 #minutes
    growth_start= pump_data[len(pump_data)-2][0] + (time_buffer/60) # hours
    growth_end = pump_data[len(pump_data)-1][0]

    ## Gets part of data that is part of the growth curve
    segmented_OD = ODData[np.where(np.logical_and(ODData >=growth_start, ODData<=growth_end))[0]]

    ## Averages the OD data from start of curve and end of curve
    averaged_values = 10
    OD_start = np.mean(segmented_OD[0:averaged_values,1])
    OD_end = np.mean(segmented_OD[len(segmented_OD)-averaged_values :len(segmented_OD),1])

    ## Calculates growth rate
    growth_rate = ((OD_end-OD_start)/OD_start)/(growth_end - growth_start )

    return growth_rate
