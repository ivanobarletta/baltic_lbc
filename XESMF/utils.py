import esmpy 
import xarray as xr
import numpy as np
import subprocess as sub

rEarth = 6317000.0  
rEarth = 6378137.0
rEarth = 6371229.0 # Nemo3.6

def slice_dataset(ds,box):
    # box = [lon0,lon1,lat0,lat1]
    # the box is extended by 5% of original size 
    # on each side
    lon_min,lon_max,lat_min,lat_max = box
    delta_lon = box[1]-box[0]
    delta_lat = box[3]-box[2]
    lon_min -= 0.05 * delta_lon
    lon_max += 0.05 * delta_lon
    lat_min -= 0.05 * delta_lat
    lat_max += 0.05 * delta_lat
    dsOut = ds.sel(lon=slice(lon_min,lon_max),lat=slice(lat_min,lat_max))
    return dsOut

def create_curvilinear_grid_from_ds(ds):
    # the input dataset must be like
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


    xSize = ds.sizes["x"]
    ySize = ds.sizes["y"]

    grid = esmpy.Grid(np.array([xSize,ySize]), 
            staggerloc=[esmpy.StaggerLoc.CENTER, esmpy.StaggerLoc.CORNER],
            coord_sys=esmpy.CoordSys.SPH_DEG,pole_dim=0)
    
    gridLon = grid.get_coords(0)
    print(gridLon.shape)
    gridLat = grid.get_coords(1)
    gridLonCorner = grid.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    gridLatCorner = grid.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)

    print(type(gridLon))

    gridLon[:] = ds.lon.transpose("x","y").data
    gridLat[:] = ds.lat.transpose("x","y").data
    gridLonCorner[:] = ds.lon_b.transpose("x_b","y_b").data
    gridLatCorner[:] = ds.lat_b.transpose("x_b","y_b").data
    
    return grid

# ----------- 
def create_test_grid():
    grid = esmpy.Grid(np.array([3,4]), 
                        staggerloc=[esmpy.StaggerLoc.CENTER, esmpy.StaggerLoc.CORNER],
                        coord_sys=esmpy.CoordSys.SPH_DEG,
                        num_peri_dims=1, periodic_dim=0, pole_dim=1)


    gridLon = grid.get_coords(0)
    gridLat = grid.get_coords(1)
    gridLonCorner = grid.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    gridLatCorner = grid.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)

    lon = np.linspace(-120,120,3)
    lat = np.linspace(-67.5, 67.5,4)
    lon_corner = np.arange(-180,180,120)
    lat_corner = np.linspace(-90, 90, 5)

    lonm, latm = np.meshgrid(lon, lat, indexing='ij')
    lonm_corner, latm_corner = np.meshgrid(lon_corner, lat_corner, indexing='ij')

    gridLon[:] = lonm
    gridLat[:] = latm
    gridLonCorner[:] = lonm_corner
    gridLatCorner[:] = latm_corner

    return grid

def create_test_grid2():
    nx = 6
    ny = 10
    xmin = 0
    xmax = 6
    ymin = 30
    ymax = 40
    grid = esmpy.Grid(np.array([nx,ny]), 
                        staggerloc=[esmpy.StaggerLoc.CENTER, esmpy.StaggerLoc.CORNER],
                        coord_sys=esmpy.CoordSys.SPH_DEG,
                        num_peri_dims=0, periodic_dim=0, pole_dim=1)

    gridLon = grid.get_coords(0)
    gridLat = grid.get_coords(1)
    gridLonCorner = grid.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    gridLatCorner = grid.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)

    lon_corner = np.linspace(xmin,xmax,nx+1)
    lat_corner = np.linspace(ymin,ymax,ny+1)
    lon = 0.5*(lon_corner[1:]+lon_corner[:-1])
    lat = 0.5*(lat_corner[1:]+lat_corner[:-1])

    lonm, latm = np.meshgrid(lon, lat, indexing='ij')
    lonm_corner, latm_corner = np.meshgrid(lon_corner, lat_corner, indexing='ij')

    gridLon[:] = lonm
    gridLat[:] = latm
    gridLonCorner[:] = lonm_corner
    gridLatCorner[:] = latm_corner

    return grid

def interp_vertical(ds,targetLevelsPath):

    # read file with z1,z2,z3...
    with open(targetLevelsPath) as f:
        targetLevels = np.asarray([float(i) for i in f.read().split(",")])

    kw={"fill_value": "extrapolate"}
    # do vertical interpolation
    dsOut = ds.interp(depth=targetLevels,kwargs=kw)

    return dsOut
    
def project_velocity(ds,dsAngles):
    # project geographical velocity onto
    # NEMO (i,j) coordinates system using
    # sin,cos from provided file

    gcost = dsAngles.gcost
    gsint = dsAngles.gsint
    # eliminate time_counter and create dimensions
    # (time_counter,y,x) --> (time,depth,y,x)
    nz = ds.sizes["depth"]
    gsint2 = gsint.isel(time_counter=0).expand_dims(dim={"depth":nz},axis=0).expand_dims(dim={"time":1},axis=0)
    gcost2 = gcost.isel(time_counter=0).expand_dims(dim={"depth":nz},axis=0).expand_dims(dim={"time":1},axis=0)

    #gsint2 = gsint.isel(time_counter=0).expand_dims(dim={"depth":nz},axis=0).expand_dims(dim={"time":1},axis=0)
    #gcost2 = gcost.isel(time_counter=0).expand_dims(dim={"depth":nz},axis=0).expand_dims(dim={"time":1},axis=0)

    #print("gsint2.shape", gsint2.data.shape)
    #print("gcost2.shape", gcost2.data.shape)

    uvel = ds["uo"]
    vvel = ds["vo"]
    # do the projection
    # | u_r |   |  gcost  -gsint |   |  u_g |
    # |     | = |                | * |      |
    # | v_r |   |  gsint   gcost |   |  v_g |
    #    

    #print("uvel.shape", uvel.data.shape)
    #print("vvel.shape", vvel.data.shape)

    ds["uo"].data =  gcost2.data * uvel.data - gsint2.data * vvel.data
    ds["vo"].data =  gsint2.data * uvel.data + gcost2.data * vvel.data
    
    return ds

def rename_dataset(ds,date_str="19000101"):
    #ds = ds.drop("time")
    #ds = ds.drop("time_counter")

    var_dictionary = dict({"thetao":"votemper",
                        "so":"vosaline",
                        "uo":"vozocrtx",
                        "vo":"vomecrty",
                        "zos_detided":"sossheig",
                        "time":"time_counter",
                        "depth":"deptht"})
    
    # rename of time dimension is done later 
    dim_dictionary = dict({
                        "x":"X",
                        "y":"Y",
                        "depth":"Z"
                        })
    


    ds = ds.rename_vars(var_dictionary)
    ds = ds.rename_dims(dim_dictionary)
    ds["votemper"].attrs["coordinates"] = "nav_lon nav_lat"
    ds["vosaline"].attrs["coordinates"] = "nav_lon nav_lat"
    ds["vozocrtx"].attrs["coordinates"] = "nav_lon nav_lat"
    ds["vomecrty"].attrs["coordinates"] = "nav_lon nav_lat"
    ds["sossheig"].attrs["coordinates"] = "nav_lon nav_lat"

    dsName = "trename_%s.nc" % date_str
    # this mess is to rename time dimension
    ds.to_netcdf(dsName)
    # renaming of UNLIMITED dimension does not work well in xarray
    # if I do ds.rename_dims({"time":"T"})
    # I get, when I do ncdump -h
    #   
    #    time = UNLIMITED ; // (0 currently)
    #    T = 1 ;
    # 
    #  the dimension T is just added
    #
    # Instead I want just to have 
    #    T = UNLIMITED ; // (1 currently)
    #
    # I do this with NCO

    result = sub.run(["ncrename","-O","-d","time,T",dsName,dsName], capture_output=True, text=True)
    if result.returncode != 0:
        raise("Problem with ncrename!!")

    ds.close()
    ds = xr.open_dataset(dsName)
    

    return ds

def reshape_dataset(ds):
    ds = ds.transpose("time","depth","x","y")
    return ds

def change_time_reference(ds):
    ds.time_counter.encoding["units"] = "hours since 1950-01-01 00:00:00"
    return ds

def add_TZXY(ds):
    # add coordinates X(X) = 1,2,3...nx
    #                 Y(Y) = 1,2,3...ny

    nt,nz,nx,ny = ds.sizes["T"],ds.sizes["Z"],ds.sizes["X"],ds.sizes["Y"]
    print (nt,nz,nx,ny)
    ds["T"] = np.arange(1,nt+1,dtype="int32")
    ds["T"].attrs = {"standard_name":"projection_t_coordinate","units":"1","axis":"T","grid_point":"T"}    
    ds["Z"] = np.arange(1,nz+1,dtype="int32")
    ds["Z"].attrs = {"standard_name":"projection_z_coordinate","units":"1","axis":"Z","grid_point":"T"}    
    ds["X"] = np.arange(1,nx+1,dtype="int32")     
    ds["X"].attrs = {"standard_name":"projection_x_coordinate","units":"1","axis":"X","grid_point":"T"}
    ds["Y"] = np.arange(1,ny+1,dtype="int32")
    ds["Y"].attrs = {"standard_name":"projection_y_coordinate","units":"1","axis":"Y","grid_point":"T"}

    return ds
    
def change_type(ds):
    # convert these variables to float32
    vars = ["votemper","vosaline","vozocrtx","vomecrty","deptht"]
    for var in vars:
        ds[var] = ds[var].astype("float32")

    return ds


"""
ds = xr.open_dataset("target_hgrid.nc")

grid_ds = create_curvilinear_grid_from_ds(ds)

field = esmpy.Field(grid_ds)
field.get_area()
areas = rEarth**2 * field.data

print (areas)
print (np.sum(areas)/ 1000000)
"""



"""

lons = ds.lon.data
lats = ds.lat.data
lons_corners = ds.lon_b.data
lats_corners = ds.lat_b.data

grid = create_curvilinear_grid(lons,lats,lons_corners,lats_corners)
"""
