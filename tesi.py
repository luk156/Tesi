#! /usr/bin/python

import os
import obspy as ob
import numpy as np
import scipy
from scipy import signal,io,fftpack
from obspy import signal, sac
import datetime as dt

ora = dt.datetime(year = 2013, month = 2, day = 16, hour = 22)
data_folder = "/home/matteo/Tesi-data/glaser-data/"
file_name = "16feb21-23.SAC"
#header,data,start = read_gyro_file(ora, data_folder)
#data100 = decimate_gyro_data(data)
#generate_raw_sac(data100[1], start, file_name )

start = dt.datetime(day = 16, month = 2, year = 2013, hour = 21, minute = 0, second = 0 )
stop = dt.datetime(day = 16, month = 2, year = 2013, hour =23 , minute = 0, second = 0)

generate_sac(start, stop, data_folder, file_name)
#fsagnac=np.transpose(header)[27]

#plt.figure()
#plt.title('Frequenza Sagnac')
#plt.plot(fsagnac,'r')

# plt.figure()
# plt.title('speed 100 Hertz')
# plt.plot(data100, 'g')

# plt.show()