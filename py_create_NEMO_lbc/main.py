##########################################
# this program creates Lateral Bounday Condition (LBC) 
# file for NEMO. 
# synopsis:
#   python main.py namelist.nml
# 
# Assumptions:
# 1) Input files of coarse dataset must be from CMEMS (Possibly working
#    also for general regular lon,lat datasets...)
# 2) The program requires as input an existing NEMO LBC file 
#    (namelist key: cn_coords_lbc) for the same LBC you want
#    to create.         

from io_utils import *
from interp_utils import interp_variables
from vel_utils import project_velocity
import sys

def main():
    
    print("main")

    """
    try:
        nml_path = sys.argv[1]
    except:
        print("You must provide a namelist")    
        sys.exit(0)     
    """
    
    nml_path = "test.nml"

    nml = get_nml(nml_path)

    print_namelist(nml)

    ds_coarse = get_coarse_coords(nml)

    # get target coordinates 
    lons,lats,levels = get_target_coords(nml)

    # get variables to interpolate
    varsDict = get_variables_info(nml)

    #print("vars_dict:")
    #print(varsDict)

    # do the interpolation
    dsOut = interp_variables(varsDict,lons,lats,levels,nml)

    #print("dsOut")
    #print(dsOut)

    renameDict = create_rename_dictionary(varsDict)

    dsOut = rename_ds(dsOut,renameDict)

    dsOut = project_velocity(dsOut,nml)

    write_lbc(dsOut,nml)

if __name__ == "__main__":
    main()