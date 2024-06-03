import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from create_cdo_target_grid import xT,yT,xf,yf
import xarray as xr
import numpy as np

anglesPath = "angles36_east_2.nc"
anglesPath = "2new_angles_east2.nc"
anglesPath = "neatl36_angles_east2.nc"

indexBox        = [3,13,31,41]
xTBox		= xT.isel(x=slice(indexBox[0],indexBox[1]),y=slice(indexBox[2],indexBox[3]))
yTBox		= yT.isel(x=slice(indexBox[0],indexBox[1]),y=slice(indexBox[2],indexBox[3]))

def get_poly_corners(lons,lats):
    lons = lons.data
    lats = lats.data
    idx = ((0,0),(0,-1),(-1,-1),(-1,0),(0,0))
    x = np.asarray([lons[i,j] for i,j in idx])
    y = np.asarray([lats[i,j] for i,j in idx])

    return x,y

xx,yy = get_poly_corners(xTBox,yTBox)

try:
    ds_angles = xr.open_dataset(anglesPath)
except:
    ds_angles = xr.open_dataset(anglesPath,decode_times=False)    

#gcost = ds_angles.gcost.isel(time_counter=0).data
#gsint = ds_angles.gsint.isel(time_counter=0).data

gcost = ds_angles["gcost"].isel(time_counter=0).data
gsint = ds_angles["gsint"].isel(time_counter=0).data

print (gcost**2+gsint**2)

print("gcost.shape",gcost.shape)

projj = ccrs.EuroPP()
data_crs = ccrs.PlateCarree()

#Figure 1
plt.figure()
ax = plt.axes(projection=projj)

ax.scatter(xT,yT,marker="o",color="r",s=0.2,transform=data_crs,zorder=2)
ax.scatter(xf,yf,marker="x",color="b",s=0.4,transform=data_crs,zorder=2)
ax.plot(xf,yf,color="k",linewidth=0.5,transform=data_crs,zorder=1)
ax.plot(xf.T,yf.T,color="k",linewidth=0.5,transform=data_crs,zorder=1)
# add line where transports are computed
#ax.plot(xT.isel(x=-1),yT.isel(x=-1),color="b",linewidth=0.5,transform=data_crs)

# add box 
#ax.plot(xx,yy,color="r",linewidth=1.2,transform=data_crs)


#ax.quiver(xT.data,yT.data,-gcost,-gsint,color="g",label=r"(-gcost,-gsint)($\hat{i}$)",transform=data_crs)
#ax.quiver(xT.data,yT.data, gsint,-gcost,color="b",label=r"( gsint,-gcost)($\hat{j}$)",transform=data_crs)

ax.coastlines(resolution='10m')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlabels_top = False
gl.ylabels_left = False

ax.set_extent([13,14,54.5,55.5])

plt.savefig("4curvilinear_LBC_mesh",dpi=600,bbox_inches="tight")

# -----
# Figure 2 (zoom)
# ------

plt.figure()
ax = plt.axes(projection=projj)

ax.scatter(xT,yT,marker="o",color="r",s=0.2,transform=data_crs)
ax.plot(xf,yf,color="k",linewidth=0.5,transform=data_crs)
ax.plot(xf.T,yf.T,color="k",linewidth=0.5,transform=data_crs)

#ax.quiver(xT.data,yT.data,-gcost,-gsint,color="g",label=r"(-gcost,-gsint)($\hat{i}$)",transform=data_crs)
#ax.quiver(xT.data,yT.data, gsint,-gcost,color="b",label=r"( gsint,-gcost)($\hat{j}$)",transform=data_crs)

#ax.quiver(xT.data,yT.data,  gcost, gsint,color="g",label=r"(  gcost, gsint)($\hat{i}$)",transform=data_crs)
#ax.quiver(xT.data,yT.data, -gsint, gcost,color="b",label=r"( -gsint, gcost)($\hat{j}$)",transform=data_crs)

ax.quiver(xT.data,yT.data,  gcost, -gsint,color="g",label=r"(  gcost, -gsint)($\hat{i}$)",transform=data_crs)
ax.quiver(xT.data,yT.data,  gsint,  gcost,color="b",label=r"(  gsint,  gcost)($\hat{j}$)",transform=data_crs)


ax.coastlines(resolution='10m')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlabels_top = False
gl.ylabels_left = False

ax.set_extent([13,14,54.5,54.7])

ax.legend()

plt.savefig("3curvilinear_LBC_mesh_and_normals",dpi=600,bbox_inches="tight")

#plt.show()


