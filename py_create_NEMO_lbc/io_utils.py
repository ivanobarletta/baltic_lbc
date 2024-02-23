import numpy as numpy
import xarray as xr
from os.path import isfile, exists
import f90nml
import sys

def get_nml(path):
    nml = f90nml.read(path)
    return nml

def print_namelist(nml):
    print(type(nml))
    print(nml)
    
def read_dataset(path):
    ds = xr.open_dataset(path)
    return ds

def get_coarse_coords(nml):
    # get coordinates
    print("Getting Coordinates of Coarse Dataset")
    path = nml["namsrc"]["cn_coord0"]
    if path == None:
        raise Exception("check namelist for cn_coord0")
    
    print("--> %s" % path)
    file_is_ok = False
    if isfile(path) and exists(path):
        file_is_ok = True
    print("file ok ", file_is_ok)     
    if not file_is_ok:
        raise Exception("The file %s does not exist" %path)

    ds = xr.open_dataset(path)
    print (ds.dims)
    return ds

def get_target_coords(nml):
    # get coords of target grid
    print("Getting Coordinates of Target Dataset")

    path = nml["namtgt"]["cn_coords_lbc"]
    if not isfile(path):
        print("the path of cn_coords_lbc not present")
        print("the program work only with NEMO lbc as")
        print("target grid at the moment")
        raise Exception("Problem with cn_coords_lbc")

    ds = xr.open_dataset(path)
    lons = ds.nav_lon
    lats = ds.nav_lat
    levels = nml["namtgt"]["z_tgt_levels"]
    
    return lons,lats,levels

def write_lbc(dsOut,nml):
    #Write dataset into netCDF file
    outFile = nml["namout"]["cn_fileout"]
    print()
    print("Writing lbc to file: %s" % outFile)
    dsOut.to_netcdf(outFile)

def get_variables_info(nml):
    # this function is just a wrapper
    new_dict = load_vars_from_cmems(nml)
    return new_dict

def load_vars_from_cmems(nml):
    # parse the cn_varfile list
    # to find the variables and the path of respective datasets
    vars_list = nml["namvar"]["cn_varfile"]
    vars_dict = parse_vars_list(vars_list)
    return vars_dict

def parse_vars_list(vlist):
    # creates a dictionary containing 
    # info on variables to interpolate
    nVars = len(vlist)
    varNames = []
    paths = []
    # split name:path
    for vv in vlist:
        vv2 = vv.split(":")
        varNames.append(vv2[0])
        paths.append(vv2[1])

    # convert input varnames to NEMO variable names
    newNames = [guess_new_name(old_name) for old_name in varNames ]

    # build dictionary
    new_dict = {}
    for vname,path,new_name in zip(varNames,paths,newNames):
        var_data = {"name" : vname, "path": path, "new_name":new_name}
        new_dict[vname] = var_data      

    return new_dict



def guess_new_name(inputStr):
    # associate inputVar <-> NEMO variable
    new_name = None
    if inputStr in ["votemper","temperature","Temperature","thetao"]:
        new_name = "votemper"
    if inputStr in ["vosaline","salinity","Salinity","so"]:
        new_name = "vosaline"
    if inputStr in ["uo","vozocrtx"]:
        new_name = "vozocrtx"
    if inputStr in ["vo","vomecrty"]:
        new_name = "vomecrty"
    if inputStr in ["sossheig","ssh","zos","zos_detided"]:
        new_name = "sossheig"  

    if new_name == None:
        raise Exception("Couldn't guess new name for variable: %s" % inputStr)    
    
    return new_name

def create_rename_dictionary(varsDict):
    # creates a dictionary of old names <-> new names
    # correspondence. The dictionary is used to 
    # rename variables for final dataset    
    oldNames = []
    newNames = []

    # browse dictionary and create 2 lists
    for variable_name, dic in varsDict.items():
        name = dic['name']
        new_name = dic['new_name']
        oldNames.append(name)
        newNames.append(new_name)
        #print(f"Variable name: {name}")
        #print(f"New name: {new_name}")

    """
    # add old/new dimension names to lists
    oldNames.append("lon")
    newNames.append("X")
    oldNames.append("lat")
    newNames.append("Y")
    """

    renameDict = dict(zip(oldNames,newNames))

    return renameDict

def rename_ds(ds,dictionary):
    # use the dictionary to rename variables
    try:
        dsOut = ds.rename(dictionary)
    except:
        print("There was a problem in renaming the Dataset variables")
        sys.exit(0)

    return dsOut
