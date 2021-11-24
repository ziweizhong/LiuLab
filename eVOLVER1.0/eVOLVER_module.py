import socket
import os.path
import shutil
import time
import pickle
import numpy as np
import numpy.matlib
import matplotlib
import scipy.signal
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import yaml
import hashlib
import glob

with open('config.yml') as f:
    config = yaml.load(f, yaml.FullLoader)

EVOLVER_IP = config['EVOLVER_IP']

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
SAVE_PATH = os.getcwd()
#EXP_DIR = os.path.join(SAVE_PATH, config['EXP_NAME'])
EXP_DIR = SAVE_PATH

OD_CAL_PATH = os.path.join(SAVE_PATH, 'OD_cal.txt')
TEMP_CAL_PATH = os.path.join(SAVE_PATH, 'temp_calibration.txt')
CONFIG_PATH = os.path.join(SAVE_PATH, 'config.yml')


def read_OD(vials):
    od_cal = np.genfromtxt("OD_cal.txt", delimiter=',')
    UDP_IP = EVOLVER_IP
    UDP_PORT = 5554
    MESSAGE='2125,2125,2125,2125,2125,2125,2125,2125,2125,2125,2125,2125,2125,2125,2125,2125,'
    sock = socket.socket(socket.AF_INET, # Internet
                     socket.SOCK_DGRAM) # UDP
    sock.settimeout(5)
    try:
        sock.sendto(MESSAGE, (UDP_IP, UDP_PORT))
        data, addr = sock.recvfrom(1024)
        data = data.split(',')
        for x in vials:
            try:
                #data[x] = np.real(od_cal[2,x] - ((np.log10((od_cal[1,x]-od_cal[0,x])/(float(data[x]) - od_cal[0,x])-1))/od_cal[3,x]))
                data[x] = (float(data[x])-od_cal[1,x])/od_cal[0,x]
            except ValueError:
                print "OD Read Error"
                data = 'empty'
                break
    except socket.timeout:
        print "UDP Timeout (OD)"
        data = 'empty'

    return data

def update_temp(vials,exp_name):
    temp_cal = np.genfromtxt("temp_calibration.txt", delimiter=',')
    
    MESSAGE = ""
    for x in vials:
        file_path =  "%s/temp_config/vial%d_tempconfig.txt" % (SAVE_PATH,x)
        data = np.genfromtxt(file_path, delimiter=',')
        temp_set = data[len(data)-1][1]
        temp_set = int((temp_set - temp_cal[1][x])/temp_cal[0][x])
        MESSAGE += str(temp_set)
        MESSAGE += ","

    UDP_IP = EVOLVER_IP
    UDP_PORT = 5553
    sock = socket.socket(socket.AF_INET, # Internet
                     socket.SOCK_DGRAM) # UDP
    sock.settimeout(5)
    sock.sendto(MESSAGE, (UDP_IP, UDP_PORT))
    try:
        data, addr = sock.recvfrom(1024)
        data = data.split(',')
        for x in vials:
            try:
                data[x] = (float(data[x]) * temp_cal[0][x]) + temp_cal[1][x]
            except ValueError:
                print "Temp Read Error"
                data = 'empty'
                break
    except socket.timeout:
        data = 'empty'
    return data

def fluid_command(MESSAGE, vial, elapsed_time, pump_wait, exp_name, time_on, file_write):
    UDP_IP = EVOLVER_IP
    UDP_PORT = 5552
    sock = socket.socket(socket.AF_INET, # Internet
                     socket.SOCK_DGRAM) # UDP
    sock.settimeout(5)

    file_path =  "%s/pump_log/vial%d_pump_log.txt" % (SAVE_PATH,vial)
    data = np.genfromtxt(file_path, delimiter=',')
    last_pump = data[len(data)-1][0]
    if ((elapsed_time- last_pump)*3600) >pump_wait:
        sock.sendto(MESSAGE, (UDP_IP, UDP_PORT))
        if file_write == 'y':
            text_file = open(file_path,"a+")
            text_file.write("%f,%s\n" %  (elapsed_time, time_on))
            text_file.close()

def stir_rate (MESSAGE):
    UDP_IP = EVOLVER_IP
    UDP_PORT = 5551
    sock = socket.socket(socket.AF_INET, # Internet
                     socket.SOCK_DGRAM) # UDP
    sock.settimeout(5)
    sock.sendto(MESSAGE, (UDP_IP, UDP_PORT))

def parse_data(data, elapsed_time, vials, exp_name, file_name):
    if data == 'empty':
        print "%s Data Empty! Skipping data log..." % file_name
    else:
        for x in vials:
            file_path =  "%s/%s/vial%d_%s.txt" % (SAVE_PATH,file_name,x,file_name)
            text_file = open(file_path,"a+")
            text_file.write("%f,%s\n" %  (elapsed_time, data[x]))
            text_file.close()

def initialize_exp(exp_name, vials):
    dir_path =  "%s" % (SAVE_PATH)
    exp_continue = raw_input('Continue from exisiting experiment? (y/n): ')
    if exp_continue == 'n':
        OD_read = read_OD(vials)
        exp_blank = raw_input('Callibrate vials to blank? (y/n): ')
        pump_init = raw_input('Turn off pumps for this experiment? (y/n)')
        if exp_blank == 'y':
            OD_initial = OD_read
        else:
            OD_initial = np.zeros(len(vials))

        start_time = time.time()
        if os.path.exists(dir_path):
            files = glob.glob(os.getcwd()+'/*')
            for f in files:
                if f != OD_CAL_PATH and f != TEMP_CAL_PATH and f != CONFIG_PATH:
                    try:
                        shutil.rmtree(f)
                    except:
                        try:
                            os.remove(f)
                        except:
                            print("Error removing %s!"%(f))

        script_name = os.path.join(SCRIPT_DIR,'custom_script.py')
        yml_name = os.path.join(SAVE_PATH,'config.yml')
        shutil.copy(yml_name, os.path.join(EXP_DIR,'config_yml_latest.txt'))
        shutil.copy(script_name, os.path.join(EXP_DIR,'custom_script_latest.txt'))

        os.makedirs("%s/OD" % dir_path)
        os.makedirs("%s/temp" % dir_path)
        os.makedirs("%s/temp_graph" % dir_path)
        os.makedirs("%s/pump_log" % dir_path)
        os.makedirs("%s/temp_config" % dir_path)
        os.makedirs("%s/growth_rate" % dir_path)
        os.makedirs("%s/ODset" % dir_path)

        os.makedirs("%s/drug_log" % dir_path)
        os.makedirs("%s/logs" % dir_path)
        os.makedirs("%s/PIDLog" % dir_path)
        os.makedirs("%s/offset_log" % dir_path)

        for x in vials:
            OD_path =  "%s/OD/vial%d_OD.txt" % (dir_path,x)
            text_file = open(OD_path,"w").close()

            growth_path =  "%s/growth_rate/vial%d_growth_rate.txt" % (dir_path,x)
            text_file = open(growth_path,"w")
            text_file.write("0,0\n0,0\n")
            text_file.close()

            temp_path =  "%s/temp/vial%d_temp.txt" % (dir_path,x)
            text_file = open(temp_path,"w").close()

            tempconfig_path =  "%s/temp_config/vial%d_tempconfig.txt" % (dir_path,x)
            text_file = open(tempconfig_path,"w")
            text_file.write("Experiment: %s vial %d, %r\n" % (exp_name, x, time.strftime("%c")))
            text_file.write("0,30\n")
            text_file.close()

            pump_path =  "%s/pump_log/vial%d_pump_log.txt" % (dir_path,x)
            text_file = open(pump_path,"w")
            text_file.write("Experiment: %s vial %d, %r\n" % (exp_name, x, time.strftime("%c")))
            if pump_init == 'y':
                text_file.write("99999999999999,0\n")
            else:
                text_file.write("0,0\n0,0\n")
            text_file.close()

            ODset_path =  "%s/ODset/vial%d_ODset.txt" % (dir_path,x)
            text_file = open(ODset_path,"w")
            text_file.write("Experiment: %s vial %d, %r\n" % (exp_name, x, time.strftime("%c")))
            text_file.write("0,0\n")
            text_file.close()

            drug_log_path =   "%s/drug_log/vial%d_drug_log.txt" % (dir_path,x)
            text_file = open(drug_log_path,"w")
            text_file.write("0,0\n0,0\n")
            text_file.close()

            log_path = "%s/logs/vial%d_log.txt" % (dir_path,x)
            text_file = open(log_path,"w")
            text_file.write("0,0\n0,0\n")
            text_file.close()

            log_path = "%s/PIDLog/vial%d_PIDLog.txt" % (dir_path,x)
            text_file = open(log_path,"w")
            #Log file consists of t, e(t), Int(e(t) dt, 0,t), d(e(t)/dt), Kp, Ki, Kd, output
            text_file.write("0,0,0,0,0,0,0,0,0\n")
            text_file.write("0,0,0,0,0,0,0,0,0\n")
            text_file.close()

            log_path = "%s/offset_log/vial%d_offset.txt" % (dir_path,x)
            text_file = open(log_path,"w")
            text_file.write("0,0\n0,0\n")
            text_file.close()

        pump_path =  "%s/pump_log/vial00_pump_log.txt" % (dir_path)
        text_file = open(pump_path,"w")
        text_file.write("Experiment: %s vial 00, %r\n" % (exp_name, time.strftime("%c")))
        if pump_init == 'y':
            text_file.write("99999999999999,0\n")
        else:
            text_file.write("0,0\n0,0\n")
        text_file.close()

    else:
        # Loads variables from previous state of experiment
        pickle_path =  "%s/%s.pickle" % (SAVE_PATH,exp_name)
        with open(pickle_path) as f:
            loaded_var  = pickle.load(f)
        x = loaded_var
        start_time = x[0]
        OD_initial = x[1]

        # Checks for any changes in custom_script or in config.yml and updates as necessary
        try:
            latest_script_filename = os.path.join(EXP_DIR,'custom_script_latest.txt')
            script_name = os.path.join(SCRIPT_DIR,'custom_script.py')
            latest_hash = hashlib.md5(open(latest_script_filename, 'rb').read()).hexdigest()
            script_hash = hashlib.md5(open(script_name, 'rb').read()).hexdigest()
        except:
            print('Error hashing custom_script.py!')
            latest_hash = 0
            script_hash = 1

        if latest_hash != script_hash:
            # copy current custom script and config.yml to txt file
            backup_filename = 'custom_script_{0}.txt'.format(time.strftime('%y%m%d_%H%M'))
            shutil.copy(script_name, os.path.join(EXP_DIR,backup_filename))
            shutil.copy(script_name, os.path.join(EXP_DIR,'custom_script_latest.txt'))

        # Gets the latest config.yml backup and checks if it is the same as the current config file
        try:
            latest_yml_name = os.path.join(EXP_DIR,'config_yml_latest.txt')
            yml_name = os.path.join(SAVE_PATH,'config.yml')
            latest_yml_hash = hashlib.md5(open(latest_yml_name, 'rb').read()).hexdigest()
            yml_hash = hashlib.md5(open(yml_name, 'rb').read()).hexdigest()
        except:
            print('Error hashing config.yml!')
            latest_yml_hash = 0
            yml_hash = 1

        if latest_yml_hash != yml_hash:
            # copy current config.yml to txt file
            yml_filename = 'config_yml_{0}.txt'.format(time.strftime('%y%m%d_%H%M'))
            shutil.copy(yml_name, os.path.join(EXP_DIR,yml_filename))
            shutil.copy(yml_name, os.path.join(EXP_DIR,'config_yml_latest.txt'))

    return start_time, OD_initial

def save_var(exp_name, start_time, OD_initial):
    pickle_path =  "%s/%s.pickle" % (SAVE_PATH,exp_name)
    with open(pickle_path, 'w') as f:
        pickle.dump([start_time, OD_initial], f)

def graph_data(vials, exp_name, file_name):
    for x in vials:
        file_path =  "%s/%s/vial%d_%s.txt" % (SAVE_PATH,file_name,x,file_name)
        data = np.genfromtxt(file_path, delimiter=',')
        plt.plot(data[:,0], data[:,1])
 ##      plt.ylim((-.05,.5))
        plot_path =  "%s/%s_graph/vial%d_%s.png" % (SAVE_PATH,file_name,x,file_name)
        plt.savefig(plot_path)
        plt.clf()

def calc_growth_rate(vials, exp_name,elapsed_time, OD_maintain):
    for x in vials:
        ## Grab Data and make setpoint
        file_path =  "%s/OD/vial%d_OD.txt" % (SAVE_PATH,x)
        OD_data = np.genfromtxt(file_path, delimiter=',')
        if np.shape(OD_data)[0] > 1010:
            time = OD_data[-1000:-1,0]
            raw_data = OD_data[-1000:-1,1]
            ## Smooth data
            window = 50
            mask = np.ones(window)/window
            raw_data = np.convolve(raw_data,mask, 'same')
            median_filter_range = 801;
            slope_data2 = movingslope(raw_data,10,1,10)
            slope_data2 = medfilt(slope_data2[:,0],median_filter_range);
            ## Calculate Average Growth for data
            average_data =((slope_data2[(median_filter_range/2):(np.size(slope_data2)-1-median_filter_range/2)])*3600)/OD_maintain
            ## Write Growth Rate to text
            log_path = "%s/growth_rate/vial%d_growth_rate.txt" % (SAVE_PATH,x)
            text_file = open(log_path,"a+")
            text_file.write("%d,%s\n" %  (elapsed_time, np.average(average_data)))
            text_file.close()


def movingslope(vec,supportlength,modelorder,dt):
    n = np.size(vec);

    ##now build the filter coefficients to estimate the slope
    if ((supportlength % 2) == 1):
        parity = 1 #odd parity
    else:
        parity = 0;
    s = (supportlength-parity)/2
    t = np.arange((-s+1-parity),s+1)[np.newaxis]
    t = t.transpose()
    coef = getcoef(t,supportlength,modelorder)

    ## Apply the filter to the entire vector
    f = scipy.signal.lfilter(-coef,1,vec);
    Dvec = np.zeros((np.size(vec),1))
    r = s+ np.arange(1,(n-supportlength+2))
    for x in np.arange(0,np.size(r)):
        Dvec[r[x]] = f[supportlength+x-1]

    for i in np.arange(1,s+1):
        t = np.arange(1,supportlength+1)[np.newaxis]
        t = t - i
        t = t.transpose()
        coef = getcoef(t,supportlength,modelorder);
        coef = coef[np.newaxis]
        m = vec[0:supportlength][np.newaxis]
        m = np.transpose(m)
        Dvec[i-1]= np.dot(coef,m)

        if i<(s + parity):
            t = np.arange(1,supportlength+1)[np.newaxis]
            t = t - supportlength + i -1
            t = t.transpose()
            coef = getcoef(t,supportlength,modelorder)
            coef = coef[np.newaxis]
            m = vec[n-supportlength : n][np.newaxis]
            m = np.transpose(m)
            Dvec[n-i]= np.dot(coef,m)

    Dvec = Dvec/dt
    return Dvec


def getcoef(t,supportlength,modelorder):
    a = numpy.matlib.repmat(t,1,modelorder+1)
    b = numpy.matlib.repmat(np.arange(0,modelorder+1),supportlength,1)
    c = np.power(a,b)
    pinvA = numpy.linalg.pinv(c)
    coef = pinvA[1,:]
    return coef

def medfilt (x, k):
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = np.zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
    return np.median (y, axis=1)

if __name__ == '__main__':
    temp_data = update_temp('23423423')
    print temp_data
