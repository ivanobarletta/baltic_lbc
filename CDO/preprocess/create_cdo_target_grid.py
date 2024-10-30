import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt

path = "coordinates.nc"
outPath = "target_hgrid2"

fmt1 = " %7.3f"
fmt2 = "%7.3f %7.3f %7.3f %7.3f   "

fmt1 = " %16.12f"
fmt2 = "%16.12f %16.12f %16.12f %16.12f   "

# this should print all precision
fmt1 = " %s"
fmt2 = "%s %s %s %s   "

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

print("xT.shape ",xT.shape)
print("xf.shape ",xf.shape)

gridType = "curvilinear"
xSize = xT.x.shape[0]
ySize = xT.y.shape[0]
gridSize = xSize * ySize
#xVals = []
#yVals = []
#xBounds = []
#yBounds = []

xVals_str = "xvals = "
yVals_str = "yvals = "

count = 0
#for yy in xT:
#    for xx in yy.data:
for j in range(ySize):
    for i in range(xSize-1,-1,-1):    
        count += 1
        #xVals_str += " %7.3f" % xx
        xx = xT.isel(x=i,y=j).data 
        xVals_str += fmt1 % xx                
        if count % 10 == 0:
            xVals_str += "\n        "
        if count == gridSize:
            xVals_str += "\n"

count = 0
#for yy in yT:
#    for xx in yy.data:
for j in range(ySize):
    for i in range(xSize-1,-1,-1):        
        count += 1
        #yVals_str += " %7.3f" % xx
        yy = yT.isel(x=i,y=j).data
        yVals_str += fmt1 % yy                
        if count % 10 == 0:
            yVals_str += "\n        "            
        if count == gridSize:
            yVals_str += "\n"


#xVals = [xx for yy in xT for xx in yy.data]
#yVals = [xx for yy in yT for xx in yy.data]

xBounds_str = "xbounds = "
yBounds_str = "ybounds = "

count = 0
for iy in range(ySize):
    #for ix in range(xSize):
    for ix in range(xSize-1,-1,-1):    
        x_sw = xf.isel(x=ix  ,y=iy  ).data
        x_se = xf.isel(x=ix+1,y=iy  ).data
        x_ne = xf.isel(x=ix+1,y=iy+1).data
        x_nw = xf.isel(x=ix  ,y=iy+1).data   
        y_sw = yf.isel(x=ix  ,y=iy  ).data
        y_se = yf.isel(x=ix+1,y=iy  ).data
        y_ne = yf.isel(x=ix+1,y=iy+1).data
        y_nw = yf.isel(x=ix  ,y=iy+1).data   
        #xBounds_str += "%7.3f %7.3f %7.3f %7.3f   " % (x_sw,x_se,x_ne,x_nw)
        #yBounds_str += "%7.3f %7.3f %7.3f %7.3f   " % (y_sw,y_se,y_ne,y_nw)
        xBounds_str += fmt2 % (x_sw,x_se,x_ne,x_nw)
        yBounds_str += fmt2 % (y_sw,y_se,y_ne,y_nw)
        count += 1
        if count % 3 == 0:
            xBounds_str += "\n        "
            yBounds_str += "\n        "
        if count == gridSize:
            xBounds_str += "\n"
            yBounds_str += "\n"


"""
fig = plt.figure()
plt.scatter(xT,yT,marker="o",color="b",s=0.2)
plt.plot(xf,yf,color="k")
plt.plot(xf.T,yf.T,color="k")
plt.show()

"""

# create grid file for CDO


f = open(outPath,"w")
f.write("gridtype = %s\n" % gridType)
f.write("gridsize = %d\n" % gridSize)
f.write("xsize    = %d\n" % xSize)
f.write("ysize    = %d\n" % ySize)
f.write(xVals_str)
f.write(xBounds_str)
f.write(yVals_str)
f.write(yBounds_str)
f.close()


