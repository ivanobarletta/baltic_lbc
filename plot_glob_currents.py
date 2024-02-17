import numpy as np
import matplotlib.pylab as plt
import xarray as xr 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from lbc_utils import get_nemo_transect_coords

""" 
transect_path = "neatl_lbc_files/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"
transect_name = "east_2"

nemo_currents_path = "cmems_glob0.083/202201/cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m_1703089728844.nc"

lbc_currents_path = "neatl_lbc_files/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"

baltic_currents_path = "cmems_baltic/BAL-NEMO_PHY-DailyMeans-20220115.nc"
"""

transect_path = "neatl_lbc_files/NEATL36_obcdta_east_1_20220115P01_R20220123.nc"
transect_name = "east_1"

nemo_currents_path = "cmems_glob0.083/202201/cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m_1703173871242.nc"

lbc_currents_path = "neatl_lbc_files/NEATL36_obcdta_east_1_20220115P01_R20220123.nc"

depth_level = 0

x_int, y_int = get_nemo_transect_coords(transect_path,transect_name)
     

# interpolate global dataset onto transect coordinates and get (u,v)
ds = xr.open_dataset(nemo_currents_path)
ds2 = ds.isel(time=0,depth=depth_level).interp(longitude=x_int,latitude=y_int) # cannot use cubic...
print("glob time ",ds.time)
print()
ds.close()

u_glob = ds2.uo
v_glob = ds2.vo 

ds2.close()

""" 
#interpolate baltic dataset onto transect as get (u,v)
ds = xr.open_dataset(baltic_currents_path)
print("baltic time", ds.time.isel(time=0))
ds2 = ds.isel(time=0,depth=depth_level).interp(lon=x_int,lat=y_int)
ds.close()
u_balt = ds2.uo 
v_balt = ds2.vo
ds2.close() 
"""

# get the data from NEMO lbc files

ds = xr.open_dataset(lbc_currents_path)
print("lbc time ",ds.time_counter.data)
u_lbc = ds.vozocrtx.isel(T=0,X=0,Z=depth_level)
v_lbc = ds.vomecrty.isel(T=0,X=0,Z=depth_level)
date_str = np.datetime_as_string(ds.time_counter.isel(T=0),unit="h")
ds.close()

print(u_lbc.shape)
### make plot
projj = ccrs.EuroPP()
data_crs = ccrs.PlateCarree()

#plt.figure()
ax = plt.axes(projection=projj)

#ax.scatter(x_int,y_int,transform=data_crs)
glob_quiver = ax.quiver(x_int,y_int,u_glob,v_glob,transform=data_crs,width=0.002,label="glob",color="c")
lbc_quiver  = ax.quiver(x_int,y_int,u_lbc ,v_lbc ,transform=data_crs,width=0.002,label="lbc",color="r")
#balt_quiver = ax.quiver(x_int,y_int,u_balt,v_balt,transform=data_crs,width=0.002,label="balt",color="y")

ax.coastlines(resolution='10m')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlabels_top = False
gl.ylabels_left = False

#ax.set_extent([13,14.5,54.5,55.5])
ax.set_extent([6,12,36,45])

figtitle = "Level %d - Time: %s" % (depth_level,date_str) 

plt.title(figtitle)
plt.legend()

filename = "currents_%s_level_%d_time_%s" % (transect_name,depth_level,date_str)

plt.savefig("currents/"+filename,dpi=600)

plt.show()

