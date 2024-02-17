import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 

neatlBathyPath  = "bathy_meter2.nc"
balticBathyPath = "bathy_BALTIC_ger_swe.nc"
neatlEast2LBC   = "../neatl_lbc_files/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"

dsNeatl     = xr.open_dataset(neatlBathyPath)
dsBaltic    = xr.open_dataset(balticBathyPath)
dsEast2     = xr.open_dataset(neatlEast2LBC)

bathyNeatl  = dsNeatl["Bathymetry"]
bathyBaltic = dsBaltic["deptho"]

nav_lon     = dsNeatl["nav_lon"]
nav_lat     = dsNeatl["nav_lat"]

lons        = dsBaltic["longitude"]
lats        = dsBaltic["latitude"]
lons,lats   = np.meshgrid(lons,lats)

nav_lon_east2 = dsEast2["nav_lon"]
nav_lat_east2 = dsEast2["nav_lat"]

def get_poly_corners(lons,lats):
    lons = lons.data
    lats = lats.data
    idx = ((0,0),(0,-1),(-1,-1),(-1,0),(0,0))
    x = np.asarray([lons[i,j] for i,j in idx])
    y = np.asarray([lats[i,j] for i,j in idx])

    return x,y

xx,yy = get_poly_corners(nav_lon_east2,nav_lat_east2)
print(xx)
print(yy)



# functions to replace default format (x,y) -> (x,y,value,jj_min,ji_min)
def fmt0(x,y):
    dist    = np.sqrt( (x - nav_lon)**2 + (y-nav_lat)**2).data
    jj,ji   = np.unravel_index(np.argmin(dist),dist.shape)
    #ji      = np.argmin(dist,axis=1)[0]     
    #jj      = np.argmin(dist,axis=0)[0]
    z       = bathyNeatl.data[jj,ji] 
    return 'x={x:.2f}  y={y:.2f} bathy={z:.2f} ji={ji:4d} jj={jj:4d}'.format(x=x, y=y, z=z, ji=ji, jj=jj)

def fmt1(x,y):
    dist    = np.sqrt( (x - lons)**2 + (y-lats)**2).data
    jj,ji   = np.unravel_index(np.argmin(dist),dist.shape)
    #ji      = np.argmin(dist,axis=1)[0]     
    #jj      = np.argmin(dist,axis=0)[0]
    z       = bathyBaltic.data[jj,ji] 
    return 'x={x:.2f}  y={y:.2f} bathy={z:.2f} ji={ji:4d} jj={jj:4d}'.format(x=x, y=y, z=z, ji=ji, jj=jj)


# make simple plot of bathymetry datasets 
vminn   = 0
vmaxx   = 50

fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(12,7))
dsNeatl["Bathymetry"].plot(ax=ax[0],x="nav_lon",y="nav_lat",vmin=vminn,vmax=vmaxx)
dsBaltic["deptho"].plot(ax=ax[1],vmin=vminn,vmax=vmaxx)
ax[0].set_xlim([12,15])
ax[0].set_ylim([54,56])
ax[0].format_coord = fmt0
ax[1].format_coord = fmt1
ax[0].set_title("NEATL36")
ax[1].set_title("CMEMS-Baltic")
ax[0].plot(xx,yy,color="k",zorder=1,linewidth=0.2)
ax[1].plot(xx,yy,color="k",zorder=1,linewidth=0.2)
ax[0].scatter(nav_lon_east2,nav_lat_east2,s=0.1,color="k")

figname = "bathy_neatl36_cmems_bal"
plt.savefig(figname,bbox_inches="tight",dpi=600)

# only neatl36
fig,ax = plt.subplots(1,1,sharex=True,sharey=True,figsize=(8,7))
dsNeatl["Bathymetry"].plot(ax=ax,x="nav_lon",y="nav_lat",vmin=vminn,vmax=vmaxx)
ax.set_xlim([12,15])
ax.set_ylim([54,56])
ax.set_title("NEATL36")
ax.scatter(nav_lon_east2,nav_lat_east2,s=0.1,color="k")

figname = "bathy_neatl36"
plt.savefig(figname,bbox_inches="tight",dpi=600)





plt.show()




