import xarray as xr 
import numpy as np
import glob
import os
import pandas as pd
from datetime import datetime
import re
import sys

"""
    Note 1)

        This script calculates transport across sections
        from "NEMO OBC like" files. 

        Files Necessary:

            - pathHGrid (containing horizontal scale factors)
            - pathZGrid (containing vertical scale factors)
            - pathUMask (containing mask)
        
            - rootPath (path where are files)
            - rootFileName (common root string for files) 
            
        the glob function builds the file list with wildcard    

    Note 2)    

        The indexing of points for OBC files

        The indexing of point along normal direction
        to the boundary are reversed (0 outmost, -1 is towards inside )

        This is not, instead, for coordinates, mask files (indexes
        with conventional order)

        that is why xidx = 0 for LBC files corresponds to
                    nx-idx-1 in static files


    Note 3)            

        l_rolling is for folders containing datasets coming from
        production systems with a date and production date.

            root_file_name_date_date_prod
            root_file_name_date_date_prod2
            root_file_name_date_date_prod3
            root_file_name_date_date_prod4

        In this case the file list contains duplicates that must
        be removed. Here the date with the latest production date
        is retained

        If the source folder contains duplicates l_rolling must be
        set to True
        
"""

pathZGrid = "../static_files/mesh_zgr_east2.nc"     # I take e3u_0 from here
pathHGrid = "../static_files/coordinates_east2.nc"  # I take e2u from here
pathUMask = "../static_files/mask_gridU_east2.nc"   # I take umask from here
xidx        = 0 # outmost index
outFolder   = "outputs"

# NEATL36 Siren LBC files
rootPath    = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/GLOBAL/PSY4V3R1/BC/TMP/"
rootFileName = "NEATL36_obcdta_east_2_*.nc"
outFileName = os.path.join(outFolder,"transport_NEATL36_east2_SIREN.nc")





# XESMF produced LBC files from CMEMS-BALtic
rootPath    = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/XESMF/outputs_bal"
rootFileName = "NEATL36_east2_BAL_XESMF_????????.nc"
outFileName = os.path.join(outFolder,"transport_NEATL36_east2_XESMF.nc")
l_rolling = False
varname = "vozocrtx"

# CDO produced LBC files from CMEMS-BALtic
rootPath    = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/CDO/outputs_bal"
rootFileName = "NEATL36_east2_BAL_CDO_????????.nc"
outFileName = os.path.join(outFolder,"transport_NEATL36_east2_CDO.nc")
l_rolling = False
varname = "vozocrtx"

# NEATL36 Siren LBC files (best analysis)
rootPath    = "/mnt/netapp1/Store_Puertos/Store_IBIop/OPERATIVAS/STORE/GLOBAL/GLO12V4/BC/TMP/"
rootFileName = "NEATL36_obcdta_east_2_*.nc"
outFileName = os.path.join(outFolder,"transport_NEATL36_east2_SIREN_best_analysis.nc")
l_rolling = True
varname = "uo"

# CDO produced LBC files from CMEMS-GLOBAL
rootPath    = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/CDO/outputs_glo"
rootFileName = "NEATL36_east2_CDO_GLO_????????.nc"
outFileName = os.path.join(outFolder,"transport_NEATL36_east2_CDO_GLO.nc")
varname = "vozocrtx"
l_rolling = False

def calcYFaces(pathZGrid=pathZGrid,pathHGrid=pathHGrid,pathUMask=pathUMask,xidx=xidx):
    # returns the areas in m**2 of the faces along a transect

    # pathZGrid must contain e3t_0
    # pathHGrid must contain e2u
    # pathUMask must contain mask

    dsZGrid = xr.open_dataset(pathZGrid)
    dsHGrid = xr.open_dataset(pathHGrid,decode_times=False)
    dsMask  = xr.open_dataset(pathUMask)

    nz = dsZGrid.sizes["z"]
    nx = dsZGrid.sizes["x"]

    # cells thickness (using value at T point, should be okay)
    h = dsZGrid["e3u_0"].isel(time=0,x=nx-xidx-1)
    # ds
    ds = dsHGrid["e2u"].isel(time=0,x=nx-xidx-1)
    # repeat ds along z (first I do isel and then I create z dim again)
    # xrray does not have "repeat" option like in numpy 
    dy = ds.isel(z=0).expand_dims(dim={"z":nz})
    # get u-mask
    umask = dsMask["mask"].isel(time_counter=0,x=nx-xidx-1).rename({"depthu":"z"})

    faces = h * dy * umask

    faces = faces.rename({"z":"Z","y":"Y"})

    #print(faces.shape)
    #area = np.nansum(faces)
    #print("area [m**2]: ",area)

    return faces

def loadDatasets(rootPath="rootPath",rootFileName="rootFileName"):
    print("Loading from folder:")
    print(rootPath)

    fileList = sorted(glob.glob(os.path.join(rootPath,rootFileName)))
    if len(fileList) == 0:
        raise Exception("loadDatasets: No Files Found!")

    # load mfdataset
    ds = xr.open_mfdataset(fileList,combine="nested",concat_dim="T")

    return ds

def loadDatasetsNoDuplicates(rootPath="rootPath",rootFileName="rootFileName"):
    """
        this function removes duplicates of dates

        there can be more than one file with referring to the 
        same date but with different production date.

        this function removes duplicates from the file list
        keeping the most recent production date

    """
    print("Loading from folder:")
    print(rootPath)

    fileList = sorted(glob.glob(os.path.join(rootPath,rootFileName)))
    if len(fileList) == 0:
        raise Exception("loadDatasets: No Files Found!")

    dates = []
    dates_prod = []
    # search for date and date_prod in file path string

    for file in fileList:
        match1 = re.search(r"\d{8}",file)
        match2 = re.search(r"R\d{8}",file)
        #print(file)
        #print("match1",match1)
        #print("match2",match2)
        date = datetime.strptime(match1.group(),"%Y%m%d")
        date_prod = datetime.strptime(match2.group()[1:],"%Y%m%d")
        dates.append(date)
        dates_prod.append(date_prod)

    nFiles = len(dates)
    dates = np.asarray(dates)
    dates_prod = np.asarray(dates_prod)
    
    if dates.shape != dates_prod.shape:
        raise Exception("problem with extraction of dates")

    # create Pandas dataframe with both dates
    df = pd.DataFrame(data = np.vstack((dates,dates_prod)).T, index=np.arange(nFiles),columns=["date","date_prod"])

    # dataFrame with dropped duplicates
    df2 = df.sort_values("date_prod").drop_duplicates("date",keep="last").sort_values("date")

    # subset fileList 
    fileListSubset = np.asarray(fileList)[df2.index].tolist()

    f = open("subset.dat","w")
    for file in fileListSubset:
        f.write(file+"\n")
    f.close()


    # load mfdataset
    ds = xr.open_mfdataset(fileListSubset,combine="nested",concat_dim="T")

    print(ds.compute())
    print()
    for dd in ds.compute().time_counter:
        print(dd)

    return ds


def calcXTransport(ds,facesArea,xidx=xidx,varname="uo"):

    # vozo (T,Z,X=xidx,Y) -> (T,Z,Y) 
    vozo = ds[varname].isel(X=xidx)

    # faces has shape (Z,Y)
    # the product will be (T,Z,Y)
    prod =  vozo.compute().data * facesArea.data    # is a numpy array

    transport = np.nansum(prod,axis=(1,2))

    return transport

def saveDataset(ds,transport,outFileName="outname.nc"):

    time = ds["time_counter"].compute()

    daOut = xr.DataArray(data=transport,dims=["time"],coords=dict(time=time.data))
    daOut = daOut.assign_attrs({"long_name":"net_eastward_transport","units":"m^3/s"})

    # create and save dataset
    dsOut = xr.Dataset()

    dsOut["transport"] = daOut

    dsOut.to_netcdf(outFileName)


def main(l_rolling = False):
    print("Calculating Face Areas of transect..")
    faces = calcYFaces(pathZGrid=pathZGrid,pathHGrid=pathHGrid,pathUMask=pathUMask,xidx=xidx)

    if l_rolling:
        print("Loading Dataset with No Duplicates...")
        ds = loadDatasetsNoDuplicates(rootPath=rootPath,rootFileName=rootFileName)
    else:
        print("Loading Dataset...")
        ds = loadDatasets(rootPath=rootPath,rootFileName=rootFileName)

    print("Calculating Transport..")
    transport = calcXTransport(ds,faces,xidx=xidx,varname=varname)

    saveDataset(ds,transport,outFileName=outFileName)


if __name__ == "__main__":
    main(l_rolling=l_rolling)
