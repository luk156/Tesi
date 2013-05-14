#! /usr/bin/python

import obspy as ob
import numpy as np
import scipy
from obspy import signal, sac

tr = ob.read('test_speed-ordine1.SAC')
#print(tr[0].stats)

tr1 = ob.read('test_speed-2.SAC')
#print(tr[0].stats)
tr1.detrend(type='linear')
tr1.filter('bandpass', freqmin=0.01,freqmax=0.1, corners=1, zerophase=1)

tr.plot()
tr1.plot()