import xarray as xr
from eofs.xarray import Eof
import numpy as np
import sys

if len(sys.argv) < 4:
    print("Error. Not enough arguments. Provide")
    print("1) inputh path")
    print("2) varname")
    print("3) # of EOFs to rebuild the field")
    sys.exit(1)

inputPath   = "zoneB_01dav_ssh_concat.nc"
varName     = "zos"
nEOFs       = 10

inputPath   = sys.argv[1]
varName     = sys.argv[2]
nEOFs       = int(sys.argv[3])


reconstructedPath   = "reconstructed_nEOFs_%d_%s" % (nEOFs, inputPath)

ssh     = xr.open_dataset(inputPath)[varName]

# remove time average 
ssh_anomaly = ssh - ssh.mean(dim="time")

# calculate weights
try:
    coslat = np.cos(np.deg2rad(ssh.coords["latitude"].values)).clip(0., 1.)
except:
    coslat = np.cos(np.deg2rad(ssh.coords["lat"].values)).clip(0., 1.)

weights = np.sqrt(coslat)[..., np.newaxis]

# make EOF analysis
solver = Eof(ssh, weights=weights)

# rebuild field
reconstructedField  = solver.reconstructedField(nEOFs)

# add again time average
reconstructedField += ssh.mean(dim="time")

# save to netCDF
reconstructedField.to_netcdf(reconstructedPath)

