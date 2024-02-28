import xarray as xr 
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np 

varname     = "so"

lon0,lat0   = 13.5,55

date        = datetime(2022,4,11)   

#ds_glo      = xr.open_dataset("../cmems_glob12/glo_anfc_subset.nc")
ds_bal      = xr.open_dataset("../cmems_baltic/bal_anfc_subset2.nc")

profile     = ds_bal[varname].interp(time=date,method="nearest").interp(lon=lon0,lat=lat0)

dpdz        = -profile.differentiate(coord="depth") # minus because depth grows towards bottom

label2      = r"$\frac{\partial (%s)}{\partial z}$" % varname

fig,ax = plt.subplots()
ax2 = ax.twiny()
ax.plot(profile,profile.depth,color="k",label="S");ax.scatter(profile,profile.depth,color="k")
ax2.plot(dpdz,dpdz.depth,color="b",label="dS/dz");ax2.scatter(dpdz,dpdz.depth,color="b")
ax.invert_yaxis()
#ax.legend()
#ax2.legend()
ax2.set_xlabel(label2,color="b",fontsize=15)
ax2.tick_params(axis="x",labelcolor="b")
ax.set_xlabel(varname,color="k",fontsize=15)
ax.tick_params(axis="x",labelcolor="k")
ax.set_ylabel("Depth [m]", fontsize = 15)
ax.grid()

halocline_depth = dpdz.depth[np.nanargmax(abs(dpdz))].data

ax.hlines(y=halocline_depth,color="0.3",xmin=-1,xmax=20)


#halocline_depth = ds_bal.depth[np.nanargmax(abs(dpdz.data).data)]

#print(halocline_depth)

plt.show()