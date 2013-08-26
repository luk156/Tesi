#! /usr/bin/python

import obspy as ob
from obspy import signal, sac, core
import datetime as dt


def read_sismo_stream(data, data_folder, starttime=None, endtime=None ):
	file_path = data_folder+data.strftime('%Y.%m.%d')+"-00.00.00.A295.*.SAC.resamp"
	print 'open', file_path
	stream = ob.read(file_path, starttime=starttime, endtime=endtime)
	return stream

sts2 = {
	'gain': 1500.0,
	'poles': [ ],
	'sensitivity': 1/1.5896e-6,
	'zeros': [] }

def read_sismo_stream2(data_folder, starttime, endtime ):
	k=starttime
	file_path = data_folder+starttime.strftime('%Y.%m.%d')+"-00.00.00.A295.*.SAC.resamp"
	print 'open', file_path
	stream = ob.read(file_path, starttime=starttime, endtime=endtime)
	for c in stream:
		c.data=c.data*1.5896e-6/1500.0
	stream.detrend('simple')
	#stream.detrend('linear')
	print stream
	k+=dt.timedelta(days = 1)
	while k<endtime+dt.timedelta(days = 1):
		file_path = data_folder+k.strftime('%Y.%m.%d')+"-00.00.00.A295.*.SAC.resamp"
		print 'open', file_path
		s = ob.read(file_path, starttime=starttime, endtime=endtime)
		for c in s:
			c.data=c.data*1.5896e-6/1500.0
		s.detrend('simple')
		#s.detrend('linear')
		stream += s
		stream.merge(method=0, interpolation_samples=0, fill_value=0)
		print stream[0].data.shape[0]
		print stream
		k+=dt.timedelta(days = 1)
	#stream.detrend(type='simple')
	return stream

def read_sismo_stream3(data_folder, starttime, endtime ):
	k=starttime
	file_path = data_folder+starttime.strftime('%Y.%m.%d')+"-00.00.00.A295.*.SAC.resamp"
	print 'open', file_path
	stream = ob.read(file_path, starttime=starttime, endtime=endtime)
	for c in stream:
		c.data=c.data*1.5896e-6/1500.0
		c.data[0]=c.data[1]
		c.data[-1]=c.data[-2]
	stream.detrend('simple')
	#stream.detrend('linear')
	print stream
	k+=dt.timedelta(days = 1)
	while k<endtime+dt.timedelta(days = 1):
		file_path = data_folder+k.strftime('%Y.%m.%d')+"-00.00.00.A295.*.SAC.resamp"
		print 'open', file_path
		s = ob.read(file_path, starttime=starttime, endtime=endtime)
		for c in s:
			c.data=c.data*1.5896e-6/1500.0
			c.data[0]=c.data[1]
			c.data[-1]=c.data[-2]
		s.detrend('simple')
		#s.detrend('linear')
		stream += s
		stream.merge()
		print stream[0].data.shape[0]
		print stream
		k+=dt.timedelta(days = 1)
	#stream.detrend(type='simple')
	return stream
