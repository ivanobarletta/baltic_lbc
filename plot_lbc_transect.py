import xarray as xr
import numpy as np
import matplotlib.pyplot as plt 


varname = "vosaline"
varname = "votemper"
varname = "vozocrtx"


nemo_lbc_path = "neatl/NEATL36_obcdta_east_2_20220515P01_R20220524.nc"
nemo_lbc_path = "neatl/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"

depth_list = [0.4940254, 1.541375, 2.645669, 3.819495, 5.078224, 6.440614, 
    7.92956, 9.572997, 11.405, 13.46714, 15.81007, 18.49556, 21.59882, 
    25.21141, 29.44473, 34.43415, 40.34405, 47.37369, 55.76429, 65.80727, 
    77.85385, 92.32607, 109.7293, 130.666, 155.8507, 186.1256, 222.4752, 
    266.0403, 318.1274, 380.213, 453.9377, 541.0889, 643.5668, 763.3331, 
    902.3393, 1062.44, 1245.291, 1452.251, 1684.284, 1941.893, 2225.078, 
    2533.336, 2865.703, 3220.82, 3597.032, 3992.484, 4405.224, 4833.291, 
    5274.784, 5727.917 ]   

ds = xr.open_dataset(nemo_lbc_path)

time = ds.time_counter
lons = ds.nav_lon
lats = ds.nav_lat.isel(X=0)
depth = ds.deptht

print(type(time))

date_str = np.datetime_as_string(ds.time_counter.isel(T=0),unit="h")

transect = ds[varname].isel(X=0,T=0)

print(transect.shape)

lats2,depth2 = np.meshgrid(lats,np.asarray(depth_list))

def contour_levels(variable):
    if variable == "votemper":
        clevs = np.linspace(5,15,41)
        ccmap = "RdYlBu_r"
        clabel = "T"
    elif variable == "vosaline":     
        clevs = np.linspace(7,15,33)
        ccmap = "RdYlBu_r"
        clabel = "PSU"    
    elif variable in ["vozocrtx","vomecrty"]:
        clevs = np.linspace(-0.35,0.35,21)
        ccmap = "RdBu"
        clabel = "m/s"

    return clevs,ccmap,clabel  

clevs,ccmap,clabel = contour_levels(varname)

vvmin = clevs[0]
vvmax = clevs[-1]

plt.figure()

#cf = plt.contourf(lats2,depth2,transect,levels=clevs,cmap=ccmap)
cf = plt.pcolormesh(lats2,depth2,transect,cmap=ccmap,vmin=vvmin,vmax=vvmax)

cb = plt.colorbar(cf)
cb.set_label(clabel)

plt.gca().invert_yaxis()
plt.ylim((50,0))

plt.xlabel("Lat Along LBC")
plt.ylabel("Depth [m]")

title = "%s - Time: %s" % (varname,date_str)

plt.title(title)
plt.savefig("neatl/%s_transect_%s" % (varname,date_str))

plt.show()


