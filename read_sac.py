import obspy as ob
import numpy as np
import scipy
from obspy import signal, sac

tr = ob.read('test_speed-3.SAC')
print(tr[0].stats)
