#! /usr/bin/python

import obspy as ob
import numpy as np
import scipy
from obspy import signal, sac
import matplotlib.pyplot as plt
import matplotlib.dates as dates

tr = ob.read('test_speed-ordine1.SAC')
#print(tr[0].stats)

tr1 = ob.read('test_speed-2.SAC')
#print(tr[0].stats)
tr1.detrend(type='linear')
tr1.filter('bandpass', freqmin=0.01, freqmax=0.1, corners=1, zerophase=1)

tr.plot()
tr1.plot()


sts2 = {
	'poles': [(-0.037+0.037j), (-0.037-0.037j) ],
	'sensitivity': 1500.0,
	'zeros': [0j, 0j] }

start = dt.datetime(day = 16, month = 2, year = 2013, hour = 21, minute = 20, second = 0 )
stop = dt.datetime(day = 16, month = 2, year = 2013, hour =21 , minute = 40, second = 0)
UTCstart = ob.UTCDateTime(start)
UTCstop = ob.UTCDateTime(stop)

tr.plot(starttime=UTCstart, endtime=UTCstop, method="full" , equal_scale= False)