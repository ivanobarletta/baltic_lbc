import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
from math import radians, cos, sin, asin, sqrt
#from index_extraction import createTransectLine,distance,haversine
from index_extraction import *

P1 = (13.55,54.5)
P2 = (13.55,55.45)

resolution  = 0.4* (1./36)  # it's better to have resolution higher than the native mesh 

lonsTransect,latsTransect = create_transect_line(P1=P1,P2=P2,resolution=resolution)

root = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/"

print("opening coordinates dataset")
pathCoords = root + "coordinates_east2.nc"

ds = xr.open_dataset(pathCoords,decode_times=False)

lonsF   = ds["glamf"].squeeze()
latsF   = ds["gphif"].squeeze()

print("extracting F indices")
indicesF = extract_f_indices(
                latsTransect=latsTransect,
                lonsTransect=lonsTransect,
                latsF=latsF,
                lonsF=lonsF,
                )

print("checking itegrity of F indices")
indicesOK = check_f_indices(indices=indicesF)

indicesV,indicesU = extract_uv_indices(indicesF=indicesF,verbose=True)

# create DataArrays with extracted F indices
tgt_x_f     = xr.DataArray(np.array([idx[1] for idx in indicesF ]), dims="f-points")
tgt_y_f     = xr.DataArray(np.array([idx[0] for idx in indicesF ]), dims="f-points")
glamf_sel   = ds["glamf"].isel(x=tgt_x_f,y=tgt_y_f).squeeze()
gphif_sel   = ds["gphif"].isel(x=tgt_x_f,y=tgt_y_f).squeeze()

# use extracted indices for U,V
tgt_y_v     = xr.DataArray(np.array([idx[0] for idx in indicesV]),dims="v-points")
tgt_x_v     = xr.DataArray(np.array([idx[1] for idx in indicesV]),dims="v-points")
glamv_sel   = ds["glamv"].isel(x=tgt_x_v,y=tgt_y_v).squeeze()
gphiv_sel   = ds["gphiv"].isel(x=tgt_x_v,y=tgt_y_v).squeeze()

tgt_y_u     = xr.DataArray(np.array([idx[0] for idx in indicesU]),dims="u-points")
tgt_x_u     = xr.DataArray(np.array([idx[1] for idx in indicesU]),dims="u-points")
glamu_sel   = ds["glamu"].isel(x=tgt_x_u,y=tgt_y_u).squeeze()
gphiu_sel   = ds["gphiu"].isel(x=tgt_x_u,y=tgt_y_u).squeeze()


"""
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.scatter(ds["glamt"],ds["gphit"],transform=ccrs.PlateCarree(),color="k",s=0.2)    # T-points
ax.coastlines(resolution="10m")
ax.scatter(ds["glamu"],ds["gphiu"],transform=ccrs.PlateCarree(),color="g",s=0.2)    # u-points
ax.scatter(ds["glamv"],ds["gphiv"],transform=ccrs.PlateCarree(),color="r",s=0.2)    # v-points
ax.scatter(ds["glamf"],ds["gphif"],transform=ccrs.PlateCarree(),color="c",s=0.4,marker="x")    # f-points
ax.plot(lonsTransect,latsTransect,color="k",transform=ccrs.PlateCarree())
# plot all cells edges
ax.plot(ds["glamf"].squeeze(),ds["gphif"].squeeze(),transform=ccrs.PlateCarree(),color="k",linewidth=0.5)    
ax.plot(ds["glamf"].squeeze().T,ds["gphif"].squeeze().T,transform=ccrs.PlateCarree(),color="k",linewidth=0.5)

ax.plot(glamf_sel,gphif_sel,color="b",transform=ccrs.PlateCarree())
numbers = np.arange(0,len(glamf_sel))+1
ax.scatter(glamf_sel,gphif_sel,color="b",transform=ccrs.PlateCarree(),s=5)
for glam,gphi,num in zip(glamf_sel,gphif_sel,numbers):
    ax.text(glam,gphi,num,color="r",transform=ccrs.PlateCarree())

ax.scatter(glamv_sel,gphiv_sel,s=15,color="g")
ax.scatter(glamu_sel,gphiu_sel,s=15,color="r")    
"""
