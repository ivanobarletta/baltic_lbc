import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt

path = "coordinates.nc"
outPath = "target_hgrid.nc"

try:
    ds = xr.open_dataset(path)
except:
    ds = xr.open_dataset(path,decode_times=False)    

ds = ds.isel(time=0,z=0)

NxDomain = ds.dims['x']
NyDomain = ds.dims['y']

NRelaxPoints = 15

xIdx1 = NxDomain-NRelaxPoints-1
xIdx2 = NxDomain-1-1

yIdx1 = 1476-1
yIdx2 = 1548-1

print ("xIdx1,xIdx2",xIdx1,xIdx2)
print ("diff", xIdx2-xIdx1)
print ("yIdx1,yIdx2",yIdx1,yIdx2)
print ("diff", yIdx2-yIdx1)

# T - points
xT = ds.glamt.isel(x=slice(xIdx1,xIdx2+1),y=slice(yIdx1,yIdx2+1))
yT = ds.gphit.isel(x=slice(xIdx1,xIdx2+1),y=slice(yIdx1,yIdx2+1))

# f - points
xf = ds.glamf.isel(x=slice(xIdx1-1,xIdx2+1),y=slice(yIdx1-1,yIdx2+1))
yf = ds.gphif.isel(x=slice(xIdx1-1,xIdx2+1),y=slice(yIdx1-1,yIdx2+1))

e1t = ds.e1t.isel(x=slice(xIdx1,xIdx2+1),y=slice(yIdx1,yIdx2+1))
e2t = ds.e2t.isel(x=slice(xIdx1,xIdx2+1),y=slice(yIdx1,yIdx2+1))




print("xT.shape ",xT.shape)
print("xf.shape ",xf.shape)

# reverse the arrays along x (NEMO east and North Open Boudndary 
# numbered from outside to inside )
xT.data = xT.data[:,::-1]
yT.data = yT.data[:,::-1]
xf.data = xf.data[:,::-1]
yf.data = yf.data[:,::-1]


ds = xr.Dataset(
        coords={
            'lon': (['y', 'x'], xT.data, {'standard_name': 'longitude'}),
            'lat': (['y', 'x'], yT.data, {'standard_name': 'latitude'}),
            'lon_b': (['y_b', 'x_b'], xf.data),
            'lat_b': (['y_b', 'x_b'], yf.data),
        }
            )


print(ds)
ds.to_netcdf(outPath)








#out_grid.to_netcdf(outPath)