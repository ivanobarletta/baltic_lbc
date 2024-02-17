import numpy as np 
import xarray as xr 
import matplotlib.pyplot as plt 
from math import sin as SIN, cos as COS, tan as TAN, sqrt as SQRT

coordPath = "coordinates.nc"

makePlot    = False

rEarth      = 6371000.0
# indexes
iT,jT       = 1000,1500
iVn,jVn     = 1000,1500
iVs,jVs     = 1000,1499

dsCoords    = xr.open_dataset(coordPath,decode_times=False)
dsCoords    = dsCoords.isel(time=0,z=0)

# coordinates of T point
zLamT   = dsCoords["glamt"].isel(x=iT,y=jT)
zPhiT   = dsCoords["gphit"].isel(x=iT,y=jT)
# coordinates of V point (north)
zLamVn  = dsCoords["glamv"].isel(x=iVn,y=jVn)
zPhiVn  = dsCoords["gphiv"].isel(x=iVn,y=jVn)
# coordinates of V point (north)
zLamVs  = dsCoords["glamv"].isel(x=iVs,y=jVs)
zPhiVs  = dsCoords["gphiv"].isel(x=iVs,y=jVs)
dsCoords.close()

dy = rEarth * np.deg2rad((zPhiVn-zPhiVs)) 
dy2 = dsCoords["e2t"].isel(x=iT,y=jT)
print(dy)
print(dy2.data)

if makePlot:
    plt.figure()
    plt.scatter(zLamT,zPhiT,marker="o",color="k",label="T")
    plt.scatter(zLamVn,zPhiVn,marker="v",color="g",label="north")
    plt.scatter(zLamVs,zPhiVs,marker="v",color="y",label="south")
    plt.legend()
    plt.xlim((zLamT-0.5,zLamT+0.5))
    plt.show()

# calculate gcost, gsint like in geo2ocean.F90
    
rad = np.pi / 180.0
rpi = np.pi 

# direction of NP in stereographic projection
zxnpt = 0. - 2. * COS( rad*zLamT ) * TAN( rpi/4. - rad*zPhiT/2. )
zynpt = 0. - 2. * SIN( rad*zLamT ) * TAN( rpi/4. - rad*zPhiT/2. )
znnpt = zxnpt*zxnpt + zynpt*zynpt
print ("zxnpt,zynpt :",zxnpt,zynpt)

# j direction segment ( V points north and South of T point) in stereographic proj
zxvvt =  2. * COS( rad*zLamVn ) * TAN( rpi/4. - rad*zPhiVn/2. ) -  2. * COS( rad*zLamVs ) * TAN( rpi/4. - rad*zPhiVs/2. )
zyvvt =  2. * SIN( rad*zLamVn ) * TAN( rpi/4. - rad*zPhiVn/2. ) -  2. * SIN( rad*zLamVs ) * TAN( rpi/4. - rad*zPhiVs/2. )
znvvt = SQRT( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
print("zxvvt,zyvvt : ",zxvvt,zyvvt)

# calculate cosine from scalar product
gcost = (zxnpt * zxvvt + zynpt * zyvvt) / znvvt
# calculate sine from vector product 
gsint = (zxnpt * zyvvt - zynpt * zxvvt) / znvvt

print (gcost,gsint)
print(gcost**2+gsint**2)



