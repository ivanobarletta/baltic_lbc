import numpy as np
import xarray as xr 
from os.path import join , isfile
from load_mesh_file import load_nemo_mesh_file

"""
compute integral of

    <z> =  \frac{\int z dv}{\int dz} 

over a region. The value of the integral depends on the region selected
"""

xidx1,xidx2 = None,None
yidx1,yidx2 = None,None

xidx1,xidx2 = 97,None
yidx1,yidx2 = 125,365

pathMesh = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/static/meshmask_SIREN/"

in_msh = 3  # mesh split in 1,2 or 3 files

ds = load_nemo_mesh_file(pathMesh=pathMesh,in_msh=in_msh)

e1t = ds["e1t"]     #(Y,X)
e2t = ds["e2t"]     #(Y,X)
e3t = ds["e3t_0"]   #(Z,Y,X)
tmask = ds["tmask"] #(Z,Y,X)
gdept_0 = ds["gdept_0"] #(Z,Y,X)

# assign coordinates to all arrays
xcoord = np.arange(e1t.sizes["X"]).astype(np.int32)+1
ycoord = np.arange(e1t.sizes["Y"]).astype(np.int32)+1
zcoord = np.arange(e3t.sizes["Z"]).astype(np.int32)+1

e1t = e1t.assign_coords({"X":xcoord,"Y":ycoord,
                    "glamt":ds["glamt"],"gphit":ds["gphit"]})

e2t = e2t.assign_coords({"X":xcoord,"Y":ycoord,
                    "glamt":ds["glamt"],"gphit":ds["gphit"]})

e3t = e3t.assign_coords({"X":xcoord,"Y":ycoord,"Z":zcoord,
                    "glamt":ds["glamt"],"gphit":ds["gphit"],"gdept_1d":ds["gdept_1d"]})

gdept_0 = gdept_0.assign_coords({"X":xcoord,"Y":ycoord,"Z":zcoord,
                    "glamt":ds["glamt"],"gphit":ds["gphit"],"gdept_1d":ds["gdept_1d"]})

tmask   = tmask.assign_coords({"X":xcoord,"Y":ycoord,"Z":zcoord,
                    "glamt":ds["glamt"],"gphit":ds["gphit"],"gdept_1d":ds["gdept_1d"]})

# compute volume
volT = e3t*e1t*e2t  # I compute in this order to get (Z,Y,X)
# mask land values
volT = volT.where(tmask==1)

volT.name = "volumeT_0"
volT = volT.assign_attrs({"units":"m^3"})

# select subArea
gdept_0 = gdept_0.isel(X=slice(xidx1,xidx2),Y=slice(yidx1,yidx2))
volT    = volT.isel(X=slice(xidx1,xidx2),Y=slice(yidx1,yidx2))

volTot = volT.sum(dim=["Z","Y","X"]).item() # [m^3]

averageDepth = (gdept_0 * volT).sum(dim=["Z","Y","X"]).item() / volTot  #[m]

daOut = xr.DataArray(data = averageDepth)

daOut.name = "averageDepth"

newAttrs = {
        "xidx1":xidx1,
        "xidx2":xidx2,
        "yidx1":yidx1,
        "yidx2":yidx2,
        "sizeX":volT.sizes["X"],
        "sizeY":volT.sizes["Y"],
        "sizeZ":volT.sizes["Z"],
        "units":"[m]"
        }

daOut = daOut.assign_attrs(newAttrs)

