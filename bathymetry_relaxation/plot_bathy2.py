import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 
import sys

neatlBathyPath  = "bathy_meter2.nc"
balticBathyPath = "bathy_BALTIC_ger_swe.nc"
neatlEast2LBC   = "../neatl_lbc_files/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"
cmemsBathyOnNeatlPath   = "CMEMSBathyOnNeatl36.nc"
relaxedBathyPath= "bathy_meter_relax2Baltic.nc"


dsNeatl     = xr.open_dataset(neatlBathyPath)
dsBaltic    = xr.open_dataset(balticBathyPath)
dsEast2     = xr.open_dataset(neatlEast2LBC)
dsRelaxed   = xr.open_dataset(relaxedBathyPath)
dsBalticOnNeatl = xr.open_dataset(cmemsBathyOnNeatlPath)

bathyNeatl  = dsNeatl["Bathymetry"]
bathyBaltic = dsBaltic["deptho"]
bathyRelaxed= dsRelaxed["Bathymetry"]
bathyBalticOnNeatl = dsBalticOnNeatl["Bathymetry"]

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

def fmt2(x,y):
    dist    = np.sqrt( (x - nav_lon)**2 + (y-nav_lat)**2).data
    jj,ji   = np.unravel_index(np.argmin(dist),dist.shape)
    #ji      = np.argmin(dist,axis=1)[0]     
    #jj      = np.argmin(dist,axis=0)[0]
    z       = bathyBalticOnNeatl.data[jj,ji] 
    return 'x={x:.2f}  y={y:.2f} bathy={z:.2f} ji={ji:4d} jj={jj:4d}'.format(x=x, y=y, z=z, ji=ji, jj=jj)

def fmt3(x,y):
    dist    = np.sqrt( (x - nav_lon)**2 + (y-nav_lat)**2).data
    jj,ji   = np.unravel_index(np.argmin(dist),dist.shape)
    #ji      = np.argmin(dist,axis=1)[0]     
    #jj      = np.argmin(dist,axis=0)[0]
    z       = bathyRelaxed.data[jj,ji] 
    return 'x={x:.2f}  y={y:.2f} bathy={z:.2f} ji={ji:4d} jj={jj:4d}'.format(x=x, y=y, z=z, ji=ji, jj=jj)


# make simple plot of bathymetry datasets 
vminn   = 0
vmaxx   = 50

fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(15,7))
cf = dsNeatl["Bathymetry"].plot(ax=ax[0],x="nav_lon",y="nav_lat",vmin=vminn,vmax=vmaxx,add_labels=False,add_colorbar=False)
dsBaltic["deptho"].plot(ax=ax[1],vmin=vminn,vmax=vmaxx,add_labels=False,add_colorbar=False)
dsBalticOnNeatl["Bathymetry"].plot(ax=ax[2],x="nav_lon",y="nav_lat",vmin=vminn,vmax=vmaxx,add_labels=False,add_colorbar=False)
dsRelaxed["Bathymetry"].plot(ax=ax[3],x="nav_lon",y="nav_lat",vmin=vminn,vmax=vmaxx,add_labels=False,add_colorbar=False)
ax[0].set_xlim([12,15])
ax[0].set_ylim([54,56])
ax[0].format_coord = fmt0
ax[1].format_coord = fmt1
ax[2].format_coord = fmt2
ax[3].format_coord = fmt3
ax[0].set_title("NEATL36")
ax[1].set_title("CMEMS-Baltic")
ax[2].set_title("CMEMS-Baltic On Neatl")
ax[3].set_title("relaxed")
ax[0].plot(xx,yy,color="k",zorder=1,linewidth=0.2)
ax[1].plot(xx,yy,color="k",zorder=1,linewidth=0.2)
ax[3].plot(xx,yy,color="k",zorder=1,linewidth=0.2)
ax[0].scatter(nav_lon_east2,nav_lat_east2,s=0.1,color="k")

cbar = fig.colorbar(cf, ax=ax[:4], shrink=0.3, location='bottom',label = "Depth [m]")

#ax[1].tick_params(labelbottom=False,labelleft=False)

figname = "bathy_neatl36_cmems_bal_relaxed"

plt.savefig(figname,bbox_inches="tight",dpi=600)

for xidx in range(-2,-15,-3):
    fig,ax = plt.subplots(figsize=(7,3))
    dsNeatl["Bathymetry"].isel(x=xidx).plot(ax=ax,x="nav_lat",label="neatl36")
    dsRelaxed["Bathymetry"].isel(x=xidx).plot(ax=ax,x="nav_lat",label="relaxed")
    dsBalticOnNeatl["Bathymetry"].isel(x=xidx).plot(ax=ax,x="nav_lat",label="BalticOnNeatl")
    ax.invert_yaxis()
    ax.legend()
    ax.set_title("X (lines inside domain) = %d" % xidx)
    ax.set_xlim(54.4,55.5)
    ax.set_ylim(50,0)
    ax.grid()
    figname = "bathy_transect_xidx%d" % xidx
    plt.savefig(figname,bbox_inches="tight",dpi=300)

#plt.show()




