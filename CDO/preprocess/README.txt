HOW TO PERFORM 3D VERTICAL INTERPOLATION WITH CDO

the scripts: 
       createDepth3DfromCMEMS.py 
       createDepth3DfromSIREN.py 

create netcdf files containing one 3D array representing 
the depths of the cells for the source and target mesh respectively.

these files are used in the 3D vertical interpolation with CDO.

the CDO command make vertical inrtpolation considering varying depths in the domain is  

                        tgtcoordinate       infile1.nc         infile2.nc        outfile.nc   
       cdo intlevelx3d,${TARGET_DEPTH3D} ${TEMPORARY_FILE} ${SOURCE_DEPTH3D} ${TEMPORARY_FILE2}
          (intlevel3d)      /\                                   /\
                            ||                                   ||
                   createDepth3DfromSIREN.py          createDepth3DfromCMEMS.py


to use the file "infile2.nc" in the vertical interpolation with intlevel3d (or 
intlevelx3d) you must regrid to the same horizontal target mesh of the input
file because when you do the vertical interpolation the horizontal has been
done. So you have to remap like is done below (you can change to another remap method..)

1)
    cdo remapcon,/path/to/target_horizontal_grid cmems_balmfc_depths.nc balmfc_depth3D_on_NEATL36_EAST2.nc


2) 

then you must extract this single array with ncks -C (do this for both tgtcoordinate and infile2.nc!)

    ncks -C -v depth3D balmfc_depth3D_on_NEATL36_EAST2.nc 2balmfc_depth3D_on_NEATL36_EAST2.nc


