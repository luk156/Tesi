#! /usr/bin/python

import obspy as ob
from obspy import signal, sac
import datetime as dt
import sismo

laser_data_folder = "/home/matteo/Tesi-data/glaser-data/"
sismo_data_folder = "/home/matteo/Tesi-data/sismometro/"

UTCstart = ob.UTCDateTime( dt.datetime(day = 4, month = 3, year = 2013, hour =3 , minute = 52, second = 0) )
UTCstop = ob.UTCDateTime( dt.datetime(day = 4, month = 3, year = 2013, hour =3 , minute = 59, second = 0) )

st=sismo.read_sismo_stream(UTCstart, sismo_data_folder, starttime=UTCstart, endtime=UTCstop)
st.detrend(type='linear')
st.differentiate()

tr = ob.read('G-Laser-2013_3_4-3.SAC', starttime=UTCstart, endtime=UTCstop)
tr.detrend(type='linear')
tr.filter('lowpass', freq=4, corners=2, zerophase=1)
tr.decimate(2)

x=np.linspace(0,360,100)

k=np.array([])

for i in x:
	a=st[1].data*np.sin((i+180)*np.pi/180)-st[2].data*np.cos((i+180)*np.pi/180)
	c=ob.signal.cross_correlation.xcorr(a,tr[0].data,2000,full_xcorr=True)
	k=np.append(k,c[0])
	print c[0], c[1]


