import xarray as xr 
from os.path import isfile, join 

def load_nemo_mesh_file(pathMesh="",in_msh=3):
    """
    pathMesh = path containing the NEMO mesh files generated with SIREN
    in_msh  = option of mesh file generation

    # in_msh = 1 -> only 1 file created (mesh_mask.nc)
    # in_msh = 2 -> 2 files created (mesh.nc and mask.nc)
    # in_msh = 3 -> 3 files created (mask.nc, mesh_hgr.nc, mesh_zgr.nc)

    """

    if in_msh == 1:
        # only 1 file mesh_mask.nc
        path    = join(pathMesh,"mesh_mask.nc")
        if not isfile(path):
            raise Exception("file %s not existing" % path)
        ds      = xr.open_dataset(path)
    if in_msh == 2:
        # 2 files: mask.nc and mesh.nc
        pathMask    = join(pathMesh,"mask.nc")
        pathMesh    = join(pathMesh,"mesh.nc")
        if not isfile(pathMask):
            raise Exception("file %s not existing" % pathMask)
        if not isfile(pathMesh):
            raise Exception("file %s not existing" % pathMesh)
        dsMask  = xr.open_dataset(pathMask)
        dsMesh  = xr.open_dataset(pathMesh)

        ds = xr.merge([dsMask,dsMesh])

    if in_msh == 3:
        # 3 files: mask.nc, mesh_hgr.nc, mesh_zgr.nc
        pathMask    = join(pathMesh,"mask.nc")
        pathMesh2D    = join(pathMesh,"mesh_hgr.nc")
        pathMesh3D    = join(pathMesh,"mesh_zgr.nc")
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

    return ds     