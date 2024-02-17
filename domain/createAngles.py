import numpy as np 
import xarray as xr 
import matplotlib.pyplot as plt 
#from math import sin as SIN, cos as COS, tan as TAN, sqrt as SQRT
from numpy import sin as SIN, cos as COS, tan as TAN, sqrt as SQRT
import sys

rpi = np.pi
rad = rpi / 180.0 

coordPath = "coordinates.nc"
outFile = "new_angles_NEATL36.nc"

dsCoords    = xr.open_dataset(coordPath,decode_times=False)
dsCoords    = dsCoords.isel(time=0,z=0)

# load coordinates
glamt = dsCoords["glamt"]
gphit = dsCoords["gphit"]
glamv = dsCoords["glamv"]
gphiv = dsCoords["gphiv"]
nav_lon = dsCoords["nav_lon"]
nav_lat = dsCoords["nav_lat"]
time = dsCoords["time"]

# direction of the NorthPole from T cells in stereographic projection
# I start calculations from the 2-nd cell along y (j index in NEMO)
zlam = glamt.isel(y=slice(1,None))
zphi = gphit.isel(y=slice(1,None))
zxnpt = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
zynpt = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
znnpt = zxnpt*zxnpt + zynpt*zynpt


#plt.figure()
#plt.quiver(zxnpt,zynpt)
#plt.show()


# j-direction: v-point segment direction (around t-point)
#zlam = glamv(ji,jj  )
#zphi = gphiv(ji,jj  )
#zlan = glamv(ji,jj-1)
#zphh = gphiv(ji,jj-1)

# I start calculations from the 2-nd cell along y (j index in NEMO)

zlam = glamv.isel(y=slice(1,None))  # coords of V points Up in T cell 
zphi = gphiv.isel(y=slice(1,None))
zlan = glamv.isel(y=slice(0,-1))    # coords of V points down in T cell 
zphh = gphiv.isel(y=slice(0,-1))

zxvvt =  2.*COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. ) - 2.*COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
zyvvt =  2.*SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. ) - 2.*SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )

znvvt = SQRT( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
znvvt = znvvt.clip(min=1e-14)   # same check as in NEMO

print("shapes")
print(zxnpt.shape)
print(zynpt.shape)
print(zxvvt.shape)
print(zyvvt.shape)
print(znvvt.shape)

# cosinus and sinus using scalar and vectorial products
gsint = -( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt  # vector product sign depends on order of v1,v2
gcost =  ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt


summ = gsint**2 + gcost**2
print ("summ min,max")
print (np.min(summ),np.max(summ))


# fill south-most line 
gcost2 = xr.zeros_like(glamt)
gsint2 = xr.zeros_like(glamt)

# filling 1st line along y=0
gcost2.isel(y=slice(1,None)).data[:] = gcost.data[:]
gcost2.isel(y=0).data[:] = gcost.isel(y=0).data[:]
gsint2.isel(y=slice(1,None)).data[:] = gsint.data[:]
gsint2.isel(y=0).data[:] = gsint.isel(y=0).data[:]

# save to netCDF
dsOut = xr.Dataset()
gcost2 = gcost2.assign_attrs({"coordinates":"nav_lon nav_lat"})
gsint2 = gsint2.assign_attrs({"coordinates":"nav_lon nav_lat"})
gcost2 = gcost2.expand_dims(dim={"time":1})
gsint2 = gsint2.expand_dims(dim={"time":1})

dsOut["gcost"] = gcost2
dsOut["gsint"] = gsint2
dsOut["nav_lon"] = nav_lon
dsOut["nav_lat"] = nav_lat
datime = dsOut["time"].assign_attrs({"units":"seconds since 2024-01-01 00:00:00"}) 
dsOut["time"] = datime
#dsOut = dsOut.rename({"time":"time_counter"})
print("----")
print(dsOut)
print("---")
print(dsOut["time"])

dsOut.to_netcdf(outFile)





