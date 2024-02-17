import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np 
import sys

"""
Method: I take a file coordinates.nc
        of a regional NEMO configuration as input.

        I need to create a curvilinear grid for xESMF
        library that accepts datasets like the one below.

        In coordinates.nc I have T-points at glamt,gphit and
        f-points at glamf,gphif that have the same size.
        I need f-points coordinates to have 1 column and 1
        row more. Each T-point must be surrounded by 4 f-points.
        The south and west T points have only 2 f-points around
        and I need to have 4 also for these. 
        I fill the additional column/row by extruding glamt,gphit

        Extrusion to west:
            I give same latitude and longitude shifted by -dx
        Extrusion to south
            I give same longitude and latitude shifted by -dy

        Note:
            constant dx,dy is an approximation but in my case
            the cells extruded will be involved in calculations
            that do not need accuracy      
"""


# the output dataset must be like
"""
dimensions:
    y = 73 ;
    x = 15 ;
    y_b = 74 ;
    x_b = 16 ;
variables:
    double lon(y, x) ;
            lon:_FillValue = NaN ;
            lon:standard_name = "longitude" ;
    double lat(y, x) ;
            lat:_FillValue = NaN ;
            lat:standard_name = "latitude" ;
    double lon_b(y_b, x_b) ;
            lon_b:_FillValue = NaN ;
    double lat_b(y_b, x_b) ;
            lat_b:_FillValue = NaN ;
    """
# where lon_b,lat_c are coordinates of corners


inputCoordsPath = "../domain/coordinates.nc"
outputCurvGrid  = "neatl36CurvGridxESMF.nc"

print("Reading Dataset: %s " % inputCoordsPath )

try:
    dsCoords = xr.open_dataset(inputCoordsPath)
except:
    dsCoords = xr.open_dataset(inputCoordsPath,decode_times=False)     

dsCoords = dsCoords.isel(time=0,z=0)

nx = dsCoords.sizes["x"]
ny = dsCoords.sizes["y"]

# acces to T points coordinates
xT = dsCoords["glamt"].data
yT = dsCoords["gphit"].data
# acces to f point coordinates
xF = dsCoords["glamf"].data
yF = dsCoords["gphif"].data

# create bigger coordinates for F
xF2 = np.zeros((ny+1,nx+1))
yF2 = np.zeros((ny+1,nx+1))

# for some reasons the these arrays are ordered like this
#   North
#   
#  (3,0   
#  (2,0)  
#  (1,0)  (1,1)
#  (0,0)  (0,1) (0,2)
#  south

print("Extruding coordinates")
# initialize extruded coordinates
xF2[1:,1:] = xF
yF2[1:,1:] = yF

dx = 1./36  # THIS IS AN APPROXIMATION BUT IS HARMELESS
dy = 1./36  # SINCE I'M NOT GOING TO INTERPOLATE BATHYMETRY
            # IN THESE EXTRUDED CELLS

# fill south row of xF2
xF2[0,1:] = xF[0,:]
# fill 1st column of xF2
xF2[1:,0] = xF[:,0] - dx

# fill south row of yF2
yF2[0,1:] = yF[0,:] - dy
# fill 1st column of yF2 
yF2[1:,0] = yF[:,0] 

# fill south-west corner of xF2,yF2
xF2[0,0] = xF2[1,0]
yF2[0,0] = yF2[0,1]

"""
# test plot
plt.figure()
plt.scatter(xT,yT,s=0.2,color="k")
#plt.scatter(xF,yF,s=0.2,color="g")
#plt.plot(xF,yF,color="g")
#plt.plot(xF.T,yF.T,color="g")
plt.plot(xF2,yF2,color="g")
plt.plot(xF2.T,yF2.T,color="g")
plt.show()
"""

print("Saving Dataset")
## save dataset
ds = xr.Dataset(
        coords={
            'lon': (['y', 'x'], xT.data, {'standard_name': 'longitude'}),
            'lat': (['y', 'x'], yT.data, {'standard_name': 'latitude'}),
            'lon_b': (['y_b', 'x_b'], xF2.data),
            'lat_b': (['y_b', 'x_b'], yF2.data),
        }
            )

ds.to_netcdf(outputCurvGrid)


