import xarray as xr 
import matplotlib.pyplot as plt 
import sys

if len(sys.argv) < 3:
    raise Exception("provide jpni,jpnj")
    sys.exit(1)


pathBathy = "bathy_meter.nc"

jpreci = 1  # number of exchange rows / cols
jprecj = 1

# number of MPI processes along i/j
jpni    = int(sys.argv[1])
jpnj    = int(sys.argv[2])

dsBathy = xr.open_dataset(pathBathy)
jpiglo  = dsBathy.sizes["x"]  
jpjglo  = dsBathy.sizes["y"]

jpi = (jpiglo - 2 * jpreci + (jpni - 1)) // jpni + 2 * jpreci
jpj = (jpjglo - 2 * jprecj + (jpnj - 1)) // jpnj + 2 * jprecj

print ("jpi", jpi)
print ("jpj", jpj)

