#! /usr/bin/python

import obspy as ob
from obspy import signal, sac, core
import datetime as dt


def read_sismo_stream(data, data_folder, starttime=None, endtime=None ):
	file_path = data_folder+data.strftime('%Y.%m.%d')+"-00.00.00.A295.*.SAC.resamp"
	print 'open', file_path
	stream = ob.read(file_path, starttime=starttime, endtime=endtime)
	for i in range(1,3):
		stream[i].stats.paz = ob.core.AttribDict({
			'poles': [(-0.037+0.037j), (-0.037-0.037j) ],
			'sensitivity': 1500.0,
			'zeros': [0j, 0j],
			'gain': 1.0})
	return stream

