#! /usr/bin/python

import os
import obspy as ob
import numpy as np
import scipy
from scipy import signal,io,fftpack
from obspy import signal, sac
import datetime as dt

# read_gyro_file(ora, data_folder)
# legge i file .DAT generati dal giroscopio e li inserisce all' interno di due vettori: header ad 1 Hertz e data a 5 KHertz
# ora e' una variabile datetime che indica data ed ora del file da aprire
# data_folder indicato il path in cui si trovano i file del giroscopio
# nel caso non siano presenti dati per un secondo questo risulta riempito di zeri

def read_gyro_file(ora, data_folder):
	header_lenght = 202
	channel_lenght = 5000
	N_channel = 4
	file_path = data_folder+str(ora.year)+"/"+str(ora.month)+"/"+str(ora.day)+"/data"+str(ora.hour)+".dat" #ricostruisco il path completo del file
	print 'open', file_path
	raw_data = np.fromfile(file_path, dtype=np.float32, count=-1) # carico il contenuto del file all' interno di un vettore
	N_seconds = raw_data.shape[0]/(channel_lenght*N_channel+header_lenght) # calcolo il numero di secondi campionati presenti nel file
	header = np.zeros([3600, header_lenght], dtype = np.float32) # inizializzo il vettore header
	data = np.zeros([N_channel, 3600 * channel_lenght], dtype = np.float32) # inizializzo il vettore data
	data_lenght = header_lenght + N_channel * channel_lenght # calcolo la lunghezza di un secondo di dati
	start = dt.datetime(
			year = raw_data[17],
			month = raw_data[16],
			day = raw_data[15],
			hour = raw_data[14],
			minute = raw_data[6],
			second = raw_data[5],
			)
	print 'start:', start," - stop:", start + dt.timedelta(hours = 1)
	# il punto iniziale del file ha esattamento lo stesso secondo e minuto dell'inizio della trace
	for j in range(N_seconds):
		sample = raw_data[ j * data_lenght : j * data_lenght + header_lenght ][47] # secondo che sto campionando
		header[sample][:] = raw_data[ j * data_lenght : j * data_lenght + header_lenght ]
		start_data = j * data_lenght + header_lenght
		i = 0
		for k in range(4):
			data[k][(sample) * channel_lenght : (sample+1) * channel_lenght ] = raw_data[ start_data + i * channel_lenght :  start_data + (i+1) * channel_lenght]
			i+=1
	return header, data, start


def decimate_gyro_data(data, low = 110, high = 200, corners = 1, zerophase = True, cal = 632.8e-9/1.35/2.0/np.pi ):
	sagn100 = scipy.signal.decimate( data[0][:], 50 )
	#a = 
	#b =
	#firstpass = lfilter(b, a, data)
    #filtered_data = lfilter(b, a, firstpass[::-1])[::-1]
	filtered_data = ob.signal.filter.bandpass( data[0][:], low, high, 5000, corners = corners, zerophase = zerophase)
	Y_hilbert = scipy.signal.hilbert(filtered_data)
	# filtro passa banda e trasformata di hilbert
	PHI = np.unwrap(np.angle(Y_hilbert)) # calcolo della fase
	speed = np.gradient(PHI)*5000*cal # calcolo della velocita' in rad/s
	speed100 = scipy.signal.decimate( speed, 50 )
	cw100 = scipy.signal.decimate( data[1][:], 50 )
	ccw100 = scipy.signal.decimate( data[2][:], 50 )
	data100 = np.vstack((sagn100, speed100, cw100, ccw100))
	return data100 # ritorno i dati a 100 Hertz

def generate_sac(start, stop, data_folder, file_name ,extra_points=1000):
	speed_trace = sac.SacIO() # inizializzo un oggetto sac
	start = start + dt.timedelta(seconds = 33)
	stop = stop + dt.timedelta(seconds = 33)
	ora = start - dt.timedelta(hours = 1) 
	header,data,start_data = read_gyro_file(ora, data_folder)
	diff_data_seconds = (start-start_data).seconds
	print diff_data_seconds
	ora = ora + dt.timedelta(hours = 1)
	speed100 = decimate_gyro_data(data)[1]
	print speed100.shape[0]
	speed100 = np.delete(speed100, range(0, diff_data_seconds*100) ) # rimuovi punti fino a start
	print speed100.shape[0]
	data_buffer = data[:, -extra_points:]
	while stop > start_data + dt.timedelta(hours = 1):
		header1,data1,start_data = read_gyro_file(ora, data_folder)
		data1 = np.append(data_buffer, data1, axis = 1)
		ora = ora + dt.timedelta(hours = 1)
		speed100_1 = decimate_gyro_data(data1)[1]
		speed100_1 = np.delete(speed100_1, range(0, extra_points/50) )
		speed100 = np.append(speed100, speed100_1)
		data_buffer = data1[:, -extra_points:]
	diff_data_seconds = (start_data + dt.timedelta(hours = 1) - stop).seconds
	speed100 = np.delete(speed100, range(speed100.shape[0] - diff_data_seconds*100, speed100.shape[0]) ) # rimuovi punti prima di stop
	print speed100.shape[0]
	speed_trace.fromarray(speed100, starttime=ob.UTCDateTime(start - dt.timedelta(seconds = 33))) # genero un traccia dall'array delle velocita'
	speed_trace.SetHvalue('kinst', 'G-Laser Pisa')
	speed_trace.SetHvalue('delta', 0.01)
	speed_trace.WriteSacBinary(file_name)
	return True

def generate_raw_sac(speed100, start, file_name ):
	speed_trace = sac.SacIO() # inizializzo un oggetto sac
	speed_trace.fromarray(speed100, starttime=ob.UTCDateTime(start)) # genero un traccia dall'array delle velocita'
	speed_trace.SetHvalue('kinst', 'G-Laser Pisa')
	speed_trace.SetHvalue('delta', 0.01)
	speed_trace.WriteSacBinary(file_name)	


