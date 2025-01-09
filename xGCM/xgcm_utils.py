import xarray as xr
import xgcm 
import numpy as np
from os.path import isfile, join
import xnemogcm
from glob import glob
import subprocess

rootPath    = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/"

pathCoords  = join(rootPath,"CONTROL_run_outs/ZNB_native_mesh/static/meshmask_SIREN/mesh_hgr.nc")
pathMesh3D  = join(rootPath,"CONTROL_run_outs/ZNB_native_mesh/static/meshmask_SIREN/mesh_zgr.nc")

def create_1d_grid(Npoints,right=True,periodic=False):
    # build a 1d xgcm grid object

    # o -  cell centers
    # | -  cell edges

    #  right == True
    #    1   2   3   4   5   6   7   8   9    
    #    o   o   o   o   o   o   o   o   o
    #      |   |   |   |   |   |   |   |   |
    #     1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5    

    #  right == False
    #    1   2   3   4   5   6   7   8   9    
    #    o   o   o   o   o   o   o   o   o
    #  |   |   |   |   |   |   |   |   |   
    # 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5     


    shift = -0.5 + right

    ds = xr.Dataset(
        coords = {
            "x_c" : (["x_c"], np.arange(1,Npoints+1)        ,{"axis":"X"},   ),        
            "x_e" : (["x_e"], np.arange(0.5,Npoints+1+shift),{"axis":"X","c_grid_axis_shift": shift}, ),
        }        
    )

    grid = xgcm.Grid(ds,periodic=periodic)

    return grid,ds  

def create_grid_from_static(points=["t","u","v","w"],
                            path2D="",
                            path3D=""):

    # convert all dimension names to lower case (X->x, Y->y ...)



    pass


def create_minimal_coords_ds(mesh_mask=xr.Dataset()):
    """Create a minimal set of coordinates from a mesh-mask dataset.

    This creates `"central"` and `"right"` grid points for the horizontal grid
    and `"central"` and `"left"` grid points in the vertical.

    """
    try:
        N_z = len(mesh_mask.coords["z"])
    except KeyError:
        N_z = len(mesh_mask.coords["nav_lev"])
    N_y = len(mesh_mask.coords["y"])
    N_x = len(mesh_mask.coords["x"])

    coords = {
        "z_c": (["z_c", ], np.arange(1, N_z + 1),
                {"axis": "Z"}),
        "z_l": (["z_l", ], np.arange(1, N_z + 1) - 0.5,
                {"axis": "Z", "c_grid_axis_shift": - 0.5}),
        "y_c": (["y_c", ], np.arange(1, N_y + 1),
                {"axis": "Y"}),
        "y_r": (["y_r", ], np.arange(1, N_y + 1) + 0.5,
                {"axis": "Y", "c_grid_axis_shift": 0.5}),
        "x_c": (["x_c", ], np.arange(1, N_x + 1),
                {"axis": "X"}),
        "x_r": (["x_r", ], np.arange(1, N_x + 1) + 0.5,
                {"axis": "X", "c_grid_axis_shift": 0.5})
    }

    return xr.Dataset(coords=coords)


def process_siren_mask(pathDir="",in_msh=3,verbose=False):
    # process mesh mask files obtained with SIREN (3.6)
    # version to be readable by xnemogcm.

    # probably the NEMO version (and SIREN accordingly) are
    # too old and they cannot be read by xnemogcm. I have to
    # make some adjustments to the files

    # SIREN allows to create mesh/mask files in 3 modes depending
    # on the namelist parameter in_msh:

    # in_msh = 1 -> only 1 file created (mesh_mask.nc)   
    # in_msh = 2 -> 2 files created (mesh.nc and mask.nc)
    # in_msh = 3 -> 3 files created (mask.nc, mesh_hgr.nc, mesh_zgr.nc)

    if verbose: print("opening static files (with in_msh: %d)" % in_msh)

    if in_msh == 1:
        # only 1 file mesh_mask.nc
        path    = join(pathDir,"mesh_mask.nc")
        if not isfile(path):
            raise Exception("file %s not existing" % path)
        ds      = xr.open_dataset(path)
    if in_msh == 2:
        # 2 files: mask.nc and mesh.nc
        pathMask    = join(pathDir,"mask.nc")
        pathMesh    = join(pathDir,"mesh.nc")        
        if not isfile(pathMask):
            raise Exception("file %s not existing" % pathMask)
        if not isfile(pathMesh):
            raise Exception("file %s not existing" % pathMesh)
        dsMask  = xr.open_dataset(pathMask)
        dsMesh  = xr.open_dataset(pathMesh)  

        ds = xr.merge([dsMask,dsMesh])
    if in_msh == 3:
        # 3 files: mask.nc, mesh_hgr.nc, mesh_zgr.nc
        pathMask    = join(pathDir,"mask.nc")
        pathMesh2D    = join(pathDir,"mesh_hgr.nc")
        pathMesh3D    = join(pathDir,"mesh_zgr.nc")                
        if not isfile(pathMask):
            raise Exception("file %s not existing" % pathMask)
        if not isfile(pathMesh2D):
            raise Exception("file %s not existing" % pathMesh2D)
        if not isfile(pathMesh3D):
            raise Exception("file %s not existing" % pathMesh3D)                
        
        dsMask      = xr.open_dataset(pathMask).drop_vars(["X","Y","Z"])
        dsMesh2D    = xr.open_dataset(pathMesh2D).drop_vars(["X","Y"])  
        dsMesh3D    = xr.open_dataset(pathMesh3D).drop_vars(["X","Y","Z"])

        ds = xr.merge([dsMask,dsMesh2D,dsMesh3D])

    # process ds
    if verbose: print("Processing mesh/mask dataset")

    # rename dimensions
    ds = ds.rename_dims({"X":"x","Y":"y","Z":"z"})

    time_dimensions = ["time_counter","t"]  # check existence of these dims

    time_dimension_exists = any([dim in ds.dims for dim in time_dimensions])

    # add time dimension
    if not time_dimension_exists:
        ds = ds.expand_dims(dim={"t":1},axis=0)

    # grab dims
    nt = ds.sizes["t"]
    nz = ds.sizes["z"]
    ny = ds.sizes["y"]
    nx = ds.sizes["x"]

    # check monotonicity of z,gdept_1d,gdepw_1d
    #  
    # (my version of SIREN is weird and puts 0 for all the 
    # values of gdept_1d for z>0). Same problem for gdepw_1d and
    # z coordinate. z coordinate is like 1,0,0,0 (non sense..)
    # I don't know what's wrong with this version of SIREN.
    # next lines of code are to fix the dataset.

    # check for z coordinate
    monotonic_z = (   ds["z"].isel(z=slice(1, None)).values 
                    - ds["z"].isel(z=slice(None,-1)).values).all() > 0
    
    monotonic = (  ds["gdept_1d"].isel(z=slice(1,None)).values 
                 - ds["gdept_1d"].isel(z=slice(None,-1)).values).all() > 0

    monotonic2 = (  ds["gdepw_1d"].isel(z=slice(1, None)).values 
                  - ds["gdepw_1d"].isel(z=slice(None,-1)).values).all() > 0

    # create new z coordinate (1,2,3...)
    new_z = xr.DataArray(data=np.arange(1,nz+1,dtype=np.int32),
                         dims=["z"],
                         coords={"z":np.arange(1,nz+1,dtype=np.int32)})

    # replace old 1d array with new
    if not monotonic_z:
        ds["z"] = new_z

    if not monotonic:
        # use the 3D info
        new_gdept_1d = ds["gdept_0"].max(dim=["y","x"])
        ds["gdept_1d"] = new_gdept_1d

    if not monotonic2:
        # use the 3D info
        new_gdepw_1d = ds["gdepw_0"].max(dim=["y","x"])
        ds["gdepw_1d"] = new_gdepw_1d

    # xnemogcm wants nav_lev..

    ds["nav_lev"] = ds["gdept_1d"].isel(t=0)
 
    return ds 

def add_z_metrics(dsIn=xr.Dataset(),verbose=False):
    # not yet tested!

    if verbose: print("adding vertical metrics if necessary")

    # If metrics lack of z factors, like in this case
    #{('X',): ['e1t', 'e1u', 'e1v', 'e1f'],
    # ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'],
    # ('Z',): []}

    xnemogcm_metrics = xnemogcm.metrics.get_metrics(dsIn)

    # I remove point progressively to create a list with the same 
    points      = ["t","u","v","w"]
    metrics_3d  = []
    for point in points:
        for key in dsIn.keys():    
            if "e3%s" % point in key:
                metrics_3d.append(key)

    if xnemogcm_metrics[('Z',)] == []:
        xnemogcm_metrics[('Z',)] = metrics_3d        

    if xnemogcm_metrics[('Z',)] == []:
        raise Warning("Warning! the dictionary of Z metrics is still empty!")

    return xnemogcm_metrics

def make_grid(dsIn=xr.Dataset(),periodic=False):

    grid = xgcm.Grid(dsIn,metrics=xnemogcm.get_metrics(dsIn),periodic=False)

    return grid


# warning!!! to make this to work I had to:
# 1) rename dimensions to lower-case (X->x..)
# 2) rename Z(z) variable to nav_lev (see line 160 in xnemogcm/domcfg.py...)

# dsdom = open_domain_cfg(datadir="/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/ 
#   STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/mask_testrun/create_meshmask",files=["mesh_mask_rename2.nc"])

# now, If I dot get_metrics(dsdom) I get

#{('X',): ['e1t', 'e1u', 'e1v', 'e1f'],
# ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'],
# ('Z',): []}

# so there's something going wrong with the vertical scale factors

# maybe I should use compute_missing_metrics?

# yes 

def open_dom_and_nemo(  pathMesh="",
                        pathNEMO="",
                        in_msh=3,
                        verbose=False,
                        newMeshFile=False):

    """
    # Wrapper to xnemogcm.open_nemo_and_domain_cfg. 
    #
    # I use a wrapper because I have first to 
    # create a mesh_mask.nc file compliant to xnemogcm 
    # standard
    """    
         
    dsMeshMask = process_siren_mask(pathDir=pathMesh,
                            in_msh=in_msh,
                            verbose=verbose)

    # save the mesh_mask.nc file in the same folder
    # of pathMesh
    dsMeshMask["standard_type"] = "xnemogcm_compliancy"
    newFilePath = join(pathMesh,"mesh_mask.nc")
    # if newMeshFile is True, remove it before creating again
    # (If the mesh file is already loaded by xnemogcm, I cannot 
    # overwrite it unless I don't delete it)

    if newMeshFile:
        if isfile(newFilePath):
            if verbose: print("removing old mesh file")
            subprocess.run("rm -f %s" % newFilePath,shell=True)

    dsMeshMask.to_netcdf(newFilePath)                         

    nemo_files_list = glob(pathNEMO)

    if len(nemo_files_list) == 0:
        raise Exception("Error! List of Nemo files is empty! check path")

    if verbose: print("opening domain and nemo output datasets")
    ds = xnemogcm.open_nemo_and_domain_cfg(
            domcfg_files=pathMesh,
            nemo_files=nemo_files_list
            )    

    metrics_new = add_z_metrics(ds,verbose=verbose)

    if verbose: print("creating xGCM grid")
    grid = xgcm.Grid(ds,metrics=metrics_new,periodic=False)

    return grid,ds 

def test_metric(grid=None,daIn=xr.DataArray(),axes=[]):

    # test the presence of scale factors

    print(grid.get_metric(daIn,axes=axes))


def compute_vorticity():
    pass 


def project_geo2nemo(ug,vg,gcost=xr.DataArray(),gsint=xr.DataArray()):
    # project velocity from geographical frame to nemo frame

    # t is the angle between the vector pointing from cell center to North pole
    # and the vector joining the v points of the cell.

    # NOTE! the input vectors must be on the cell centers!! (y_c,x_c)

    # do the projection
    # | u_p |   |  gcost  -gsint |   |  u_g |
    # |     | = |                | * |      |
    # | v_p |   |  gsint   gcost |   |  v_g |
    #

    try:
        gcost = gcost.rename({"x":"x_c","y":"y_c"})
        gsint = gsint.rename({"x":"x_c","y":"y_c"})
    except:
        pass 

    up = gcost * ug - gsint * vg
    vp = gsint * ug + gcost * vg

    return (up,vp)

def project_nemo2geo(up,vp,gcost=xr.DataArray(),gsint=xr.DataArray()):

    # project velocity from nemo frame to geographical frame

    # t is the angle between the vector pointing from cell center to North pole
    # and the vector joining the v points of the cell.

    # NOTE! the input vectors must be on the cell centers!! (y_c,x_c)

    # do the projection
    # | u_g |   |  gcost   gsint |   |  u_p |
    # |     | = |                | * |      |
    # | v_g |   | -gsint   gcost |   |  v_p |
    #

    # check on the names of dims/coordinates. I need coordinates to be
    # exactly the same to make the product
    try:
        gcost = gcost.rename({"x":"x_c","y":"y_c"})
        gsint = gsint.rename({"x":"x_c","y":"y_c"})
    except:
        pass 

    ug =  gcost * up + gsint * vp 
    vg = -gsint * up + gcost * vp     

    print("ug.sizes")    
    print(ug.sizes)

    ug = ug.transpose("t","z_c","y_c","x_c")
    vg = vg.transpose("t","z_c","y_c","x_c")

    print("ug.sizes")    
    print(ug.sizes)


    return (ug,vg)

def calc_density(daS,daT,verbose=False):

    """
    daS = dataArray of salinity    with xgcm conventions {t,z_c,y_c,x_c}
    daT = dataArray of temperature with xgcm conventions {t,z_c,y_c,x_c}    

    !NOTE! the computation is quite slow, so the time dimension (t)
    cannot be too large (1-5/10 max)
    """

    from seawater import dens

    # build 3D depth
    if "gdept_0" not in daS.coords:
        raise Exception("Error! the input dataArray must have gdept_0 coordinate")

    # use depth as pressure (1m \approx 1db)    
    if verbose:
        print("Building 4D depth array")
    depth4D = daS["gdept_0"].expand_dims(dim={"t":daS.sizes["t"]},axis=0)

    if verbose:
        print("computing density..")
    density = dens(s=daS,t=daT,p=depth4D)

    if verbose:
        print("creating output dA")    
    densitydA = xr.DataArray(
                            data = density,
                            dims=daS.dims,
                            coords=daS.coords
                            )    

    densitydA.name = "density"                           

    return densitydA

def calc_N2(daDens,daGrid,zCoord=xr.DataArray(),verbose=False):

    """
    compute Brunt-Vaisala**2 frequency [s**-2]

    daDens = dataArray of density with xgcm conventions {t,z_c,y_c,x_c}  
    daGrid = xgcm.grid.Grid corresponding to daDens

    NOTE! daDens time dimension possibly does not match with 
    the daGrid. This happens if you create the grid object from a
    ds1 dataSet

    grid = xgcm.Grid(ds1,metrics=metrics)

    but then you make take a variable from the dataset and slice it
    along time

    varDa = ds1[var].isel(t=slice(t1,t2))

    if the time dimension of varDa is different from the time dimension 
    of ds1, operations on the grid will fail

    for example:
        grid.diff(varDa,axis="Z") 

    will rais an error for mismatch of time dimension, even though the 
    differentiation is done along Z.    

    It is safer to iterate over the time dimension of daDens
    and, after, stack the resulting outcomes into a single dataArray
    
    """
    grav = 9.8  # m/s**2

    daListdRhoDz = []   # density vertical derivative (z_f)
    daListdRhoZf = []   # density at interfaces (z_f)

    if verbose:
        print("Looping through time dimension of daDens")
    for time in range(daDens.sizes["t"]):
        dRhodz = daGrid.diff(daDens.isel(t=time),axis="Z")
        RhoZf  = daGrid.interp(daDens.isel(t=time),axis="Z")
        daListdRhoDz.append(dRhodz)
        daListdRhoZf.append(RhoZf)

    # concatenate dataArrays
    if verbose:
        print("Concatenating DataArrays")
    dRhodz = xr.concat(daListdRhoDz,dim="t")      
    RhoZf  = xr.concat(daListdRhoZf,dim="t")

    daN2 = grav * (dRhodz / RhoZf)

    daN2.name = "N"
    daN2 = daN2.assign_attrs({"units":"1/s**2","long_name":"Brunt-Vaisala frequency (square)"})
    daN2 = daN2.assign_coords({"gdepw_1d":zCoord})

    return daN2 
