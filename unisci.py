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

def read_gyro_file(ora, data_folder, header_lenght =202, channel_array = np.array([5000,5000,5000,5000])):
	ora = dt.datetime(
			year = ora.year,
			month = ora.month,
			day = ora.day,
			hour = ora.hour,
			)
	file_path = data_folder+str(ora.year)+"/"+str(ora.month)+"/"+str(ora.day)+"/data"+str(ora.hour)+".dat" #ricostruisco il path completo del file
	print 'open', file_path
	raw_data = np.fromfile(file_path, dtype=np.float32, count=-1) # carico il contenuto del file all' interno di un vettore
	N_seconds = raw_data.shape[0]/(channel_array.sum()+header_lenght) # calcolo il numero di secondi campionati presenti nel file
	print "il file contiene ", N_seconds, " secondi" 
	header = np.zeros([3600, header_lenght], dtype = np.float32) # inizializzo il vettore header
	data=[]
	for k in range(channel_array.shape[0]):
		data.append(np.zeros(channel_array[k]*3600))
	data_lenght = header_lenght + channel_array.sum() # calcolo la lunghezza di un secondo di dati
	start_trace = dt.datetime(
			year = raw_data[10],
			month = raw_data[9],
			day = raw_data[8],
			hour = raw_data[7],
			minute = raw_data[6],
			second = raw_data[5],
			)
	if start_trace>ora:
		print "il file contiene l'inizio della traccia pertanto potrebbero esserci degli artefatti"
	for j in range(N_seconds):
		sample = raw_data[ j * data_lenght : j * data_lenght + header_lenght ][47] # secondo che sto campionando
		header[sample][:] = raw_data[ j * data_lenght : j * data_lenght + header_lenght ]
		start_data = j * data_lenght + header_lenght
		i = 0
		for k in range(4):
			data[k][(sample) * channel_array[k] : (sample+1) * channel_array[k] ] = raw_data[ start_data + i :  start_data + i + channel_array[k]]
			i+=channel_array[k]
	return header, data, ora - dt.timedelta(seconds=33)


def decimate_gyro_data(sagnac, low = 110, high = 200, corners = 1, zerophase = True, cal = 632.8e-9/1.35/2.0/np.pi ):
	filtered_data = ob.signal.filter.bandpass( sagnac, low, high, 5000, corners = corners, zerophase = zerophase)
	Y_hilbert = scipy.signal.hilbert(filtered_data)
	# filtro passa banda e trasformata di hilbert
	PHI = np.unwrap(np.angle(Y_hilbert)) # calcolo della fase
	speed = np.gradient(PHI)*5000*cal # calcolo della velocita' in rad/s
	speed100 = scipy.signal.decimate( speed, 50 )
	return speed100 # ritorno la velocita' a 100 Hertz

def generate_sac(start, stop, data_folder, destination_folder="./", file_name='default' ,extra_points=1000, channel_array = np.array([5000,5000,5000,5000]), sagnac_index = 1):
	speed_trace = sac.SacIO() # inizializzo un oggetto sac
	ora = start - dt.timedelta(seconds=33+extra_points/2500) 
	header,data,start_data = read_gyro_file(ora, data_folder, channel_array = channel_array )
	diff_data_seconds = (start-start_data).seconds
	ora = ora + dt.timedelta(hours = 1)
	speed100 = decimate_gyro_data( data[sagnac_index] )
	speed100 = np.delete( speed100, range(0, diff_data_seconds*100) ) # rimuovi punti fino a start
	data_buffer = data[sagnac_index][-extra_points:]
	while stop > start_data + dt.timedelta(hours = 1):
		header1,data1,start_data = read_gyro_file(ora, data_folder,  channel_array = channel_array )
		data1 = np.append( data_buffer, data1[sagnac_index])
		ora = ora + dt.timedelta(hours = 1)
		speed100_1 = decimate_gyro_data(data1 )
		speed100_1 = np.delete(speed100_1, range(0, extra_points/50) )
		speed100 = np.append(speed100, speed100_1)
		data_buffer = data1[-extra_points:]
	diff_data_seconds = (start_data + dt.timedelta(hours = 1, seconds=33) - stop).seconds
	speed100 = np.delete(speed100, range(speed100.shape[0] - diff_data_seconds*100, speed100.shape[0]) ) # rimuovi punti prima di stop
	speed_trace.fromarray(speed100, starttime=ob.UTCDateTime(start)) # genero una traccia dall'array delle velocita'
	speed_trace.SetHvalue('kinst', 'G-Laser Pisa')
	speed_trace.SetHvalue('delta', 0.01)
	if file_name=="default":
		file_name="G-Laser-"+str(start.year)+"_"+str(start.month)+"_"+str(start.day)+"-"+str(start.hour)+":"+str(start.minute)+"_"+str(stop.hour)+":"+str(stop.minute)+".SAC"
	print "writing ", file_name
	os.chdir(destination_folder)
	speed_trace.WriteSacBinary(file_name)
	return file_name

def generate_raw_sac(speed100, start, destination_folder="./", file_name='default' ):
	speed_trace = sac.SacIO() # inizializzo un oggetto sac
	speed_trace.fromarray(speed100, starttime=ob.UTCDateTime(start)) # genero un traccia dall'array delle velocita'
	speed_trace.SetHvalue('kinst', 'G-Laser Pisa')
	speed_trace.SetHvalue('delta', 0.01)
	if file_name=="default":
		file_name="G-Laser-"+str(start.year)+"_"+str(start.month)+"_"+str(start.day)+"-"+str(start.hour)+".SAC"
	os.chdir(destination_folder)
	speed_trace.WriteSacBinary(file_name)
	print "creato il file:", file_name


