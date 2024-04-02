import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt
import sys
from numpy.fft import fft

if len(sys.argv) < 2:
    print("Error: Provide path of .nc with model ssh")
    sys.exit(1)

path    = sys.argv[1]

ds      = xr.open_dataset(path)
ssh     = ds["ssh"].squeeze()
time    = ds["time_counter"].data
ds.close()

# compute fft
fft_ssh = fft(ssh)

N       = len(fft_ssh)
n       = np.arange(N)

# sampling rate
dt      = (time[1]-time[0]) / np.timedelta64(1,"s")
sr      = 1. / dt 

T       = N / sr
freq    = n / T 

# take 1 side of fft
n_oneside   = N // 2

f_oneside   =  freq[:n_oneside]

t_days      = 1/f_oneside / 86400.0

# make plot
plt.figure()
plt.loglog(t_days, np.abs(fft_ssh[:n_oneside])/n_oneside)
plt.xlabel("cycles per day")
plt.grid()

plt.show()
