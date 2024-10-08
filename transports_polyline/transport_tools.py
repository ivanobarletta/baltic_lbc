import numpy as np
import xarray as xr
from index_extraction import *
from glob import glob
from os.path import isfile
import sys

def calcNEMOTransport(                                  # Wrapper to the main function
                        P1 = (),                        # start point (lon,lat) tuple
                        P2 = (),                        # end point (lon,lat) tuple
                        resolution = 1./36,             # resolution (in Deg) of transect line (tip: multiply the mesh resolution by a factor 0.4 -> resolution = 0.4*mesh_resolution)
                        pathCoords2D = "",              # path of coordinates.nc # e1u,e2u,e1v,e2v
                        pathZMesh = "",                 # path of file with vertical scale factors   #e3t
                        pathU = "",                     # path with U velocity files (you can use wildcards)
                        pathV = "",                     # path with V velocity files (you can use wildcards)  
                        pathS = "",                     # path of Salinity files (you can use wildcards)
                        pathT = "",                     # path of Salinity files (you can use wildcards)
                        pathMaskU = "",                 
                        pathMaskV = "",                        
                        outFileRoot = "out_transport_%.nc", # outfile root
                        transectName = "",              # name of transect (for outFile)
                        computeS = False,               # compute salinity transport (if True you must set pathS)
                        varNameS = "so",                # name of salinity variable in fileS       
                        computeH = False,               # compute heat transport (if True you must set pathT)
                        varNameT = "thetao",            # name of temperature variable in fileT       
                        verboseLevel = 0,
                        makePlot = False                          
                    ):
    
    if verboseLevel>0:
        print("###########################################")        
        print("############ calcNEMOTransport ############")
        print("###########################################")
        print("    verboseLevel: %d" % verboseLevel)
        print("    transectName: %s" % transectName)
        print("    compute Salt Transport: %s" % computeS)
        print("    compute Heat Transport: %s" % computeH)

    rhoWater            = 1025.0    # sea water density [kg/m^3]
    cpWater             = 3996.0    # sea water specific heat [J/K/Kg]         
    volumeConversion    = 1e-6      # m^3/s -> Sv
    saltConversion      = 1e-12      
    heatConversion      = 1e-15
    tempRef             = 0         # reference temperature      

    lonsTransect,latsTransect = create_transect_line(P1=P1,P2=P2,resolution=resolution)

    if verboseLevel>1:
        print(lonsTransect)
        print(latsTransect)

    # open 2D coordinates dataset
    if not isfile(pathCoords2D):
        raise Exception("Error: coords2D file not existing: %s " % pathCoords2D)
    dsCoords = xr.open_dataset(pathCoords2D,decode_times=False)

    # open mesh Z dataset
    if not isfile(pathZMesh):
        raise Exception("Error: meshZ file not existing: %s " % pathZMesh)
    dsMeshZ = xr.open_dataset(pathZMesh)

    # open U mask dataset
    if not isfile(pathMaskU):
        raise Exception("Error: maskU file not existing")
    dsMaskU = xr.open_dataset(pathMaskU)

    # open V mask dataset
    if not isfile(pathMaskV):
        raise Exception("Error: maskV file not existing")
    dsMaskV = xr.open_dataset(pathMaskV)

    # load coordinates of mesh fpoints
    lonsF   = dsCoords["glamf"].squeeze()
    latsF   = dsCoords["gphif"].squeeze()


    """
    if makePlot:
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.scatter(dsCoords["glamt"],dsCoords["gphit"],transform=ccrs.PlateCarree(),color="k",s=0.2)    # T-points
        ax.coastlines(resolution="10m")
        ax.scatter(dsCoords["glamu"],dsCoords["gphiu"],transform=ccrs.PlateCarree(),color="g",s=0.2)    # upoints
        ax.scatter(dsCoords["glamv"],dsCoords["gphiv"],transform=ccrs.PlateCarree(),color="r",s=0.2)    # vpoints
        ax.scatter(dsCoords["glamf"],dsCoords["gphif"],transform=ccrs.PlateCarree(),color="c",s=0.4,marker="x")    # f-points
        ax.plot(lonsTransect,latsTransect,color="k",transform=ccrs.PlateCarree(),linewidth=0.7)
        ax.scatter(lonsTransect,latsTransect,color="k",transform=ccrs.PlateCarree(),marker="o")

        # plot all cells edges
        ax.plot(dsCoords["glamf"].squeeze(),dsCoords["gphif"].squeeze(),transform=ccrs.PlateCarree(),color="k",linewidth=0.5)    
        ax.plot(dsCoords["glamf"].squeeze().T,dsCoords["gphif"].squeeze().T,transform=ccrs.PlateCarree(),color="k",linewidth=0.5)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
            linewidth=1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.ylabels_left = False

        plt.savefig("polyline_%s" % transectName,dpi=600,bbox_inches="tight")
        plt.show()
    """

    # extract indices of F-points
    indicesF = extract_f_indices(
                latsTransect=latsTransect,
                lonsTransect=lonsTransect,
                latsF=latsF,
                lonsF=lonsF,
                verboseLevel=verboseLevel
                )

    if verboseLevel>0:
        print("indicesF")
        print(indicesF)

    # check integrity of F indices
    indicesOK = check_f_indices(indices=indicesF,verboseLevel=verboseLevel)

    # extract U,V indices from F indices
    indicesV,indicesU = extract_uv_indices(indicesF=indicesF,verboseLevel=verboseLevel)

    # initialize nSegments along U,V
    nSegmentsU  = len(indicesU)
    nSegmentsV  = len(indicesV)
    nSegments   = nSegmentsU + nSegmentsV

    if verboseLevel>0:
        print("indicesU")
        print(indicesU)
        print("indicesV")
        print(indicesV)

    # Do some checks on the number of segments
    transectType = "general"
    if nSegments == nSegmentsU:
        # the polyline that approximates the transect is made only of U points (dy)
        #  
        #      .----.
        # j+1  |u T |u       u*dy
        #      .----.
        # j    |u T |u       u*dy
        #      .----.
        # j-1  |u T |u       u*dy
        #      .----.
        #
        #        i
        #

        print("Purely Y type transect")
        transectType = "pureY"
    if nSegments == nSegmentsV:
        # the polyline that approximates the transect is made only of V points (dx)
        #         
        #        v    v    v    v
        #     .----.----.----.----.
        #  j  |  T |  T |  T |  T |
        #     .----.----.----.----. 
        #       i-1   i   i+1  i+2
        #
        #       v*dx + v*dx ...
        #
        print("Purely X type transect")
        transectType = "pureX"

    if makePlot:
        print("making plot")
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

        # create DataArrays with extracted F indices
        tgt_x_f     = xr.DataArray(np.array([idx[1] for idx in indicesF ]), dims="fpoints")
        tgt_y_f     = xr.DataArray(np.array([idx[0] for idx in indicesF ]), dims="fpoints")
        glamf_sel   = dsCoords["glamf"].isel(x=tgt_x_f,y=tgt_y_f).squeeze()
        gphif_sel   = dsCoords["gphif"].isel(x=tgt_x_f,y=tgt_y_f).squeeze()

        # use extracted indices for U,V
        if transectType != "pureY":
            tgt_y_v     = xr.DataArray(np.array([idx[0] for idx in indicesV]),dims="vpoints")
            tgt_x_v     = xr.DataArray(np.array([idx[1] for idx in indicesV]),dims="vpoints")
            glamv_sel   = dsCoords["glamv"].isel(x=tgt_x_v,y=tgt_y_v).squeeze()
            gphiv_sel   = dsCoords["gphiv"].isel(x=tgt_x_v,y=tgt_y_v).squeeze()

        if transectType != "pureX":
            tgt_y_u     = xr.DataArray(np.array([idx[0] for idx in indicesU]),dims="upoints")
            tgt_x_u     = xr.DataArray(np.array([idx[1] for idx in indicesU]),dims="upoints")
            glamu_sel   = dsCoords["glamu"].isel(x=tgt_x_u,y=tgt_y_u).squeeze()
            gphiu_sel   = dsCoords["gphiu"].isel(x=tgt_x_u,y=tgt_y_u).squeeze()

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.scatter(dsCoords["glamt"],dsCoords["gphit"],transform=ccrs.PlateCarree(),color="k",s=0.2,label="T")    # T-points
        ax.coastlines(resolution="10m")
        ax.scatter(dsCoords["glamu"],dsCoords["gphiu"],transform=ccrs.PlateCarree(),color="g",s=0.2,label="U")    # upoints
        ax.scatter(dsCoords["glamv"],dsCoords["gphiv"],transform=ccrs.PlateCarree(),color="r",s=0.2,label="V")    # vpoints
        ax.scatter(dsCoords["glamf"],dsCoords["gphif"],transform=ccrs.PlateCarree(),color="c",s=0.4,marker="x")    # f-points
        ax.plot(lonsTransect,latsTransect,color="k",transform=ccrs.PlateCarree(),linewidth=0.7)
        # plot all cells edges
        ax.plot(dsCoords["glamf"].squeeze(),dsCoords["gphif"].squeeze(),transform=ccrs.PlateCarree(),color="k",linewidth=0.5)    
        ax.plot(dsCoords["glamf"].squeeze().T,dsCoords["gphif"].squeeze().T,transform=ccrs.PlateCarree(),color="k",linewidth=0.5)

        ax.plot(glamf_sel,gphif_sel,color="b",transform=ccrs.PlateCarree(),linewidth=0.7)
        numbers = np.arange(0,len(glamf_sel))+1
        ax.scatter(glamf_sel,gphif_sel,color="b",transform=ccrs.PlateCarree(),s=2)
        for glam,gphi,num in zip(glamf_sel,gphif_sel,numbers):
            ax.text(glam,gphi,num,color="y",transform=ccrs.PlateCarree(),fontsize=3)

        if transectType != "pureY":
            ax.scatter(glamv_sel,gphiv_sel,s=0.5,color="g")
        if transectType != "pureX":     
            ax.scatter(glamu_sel,gphiu_sel,s=0.5,color="r")   

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
            linewidth=1, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.ylabels_left = False

        plt.legend()

        plt.savefig("polyline_%s" % transectName,dpi=600,bbox_inches="tight")
        plt.show()


    # create list of files
    listFileU = sorted(glob(pathU))
    listFileV = sorted(glob(pathV))

    # fill salinity file list (list of None if computeS == False)
    if computeS:
        listFileS = sorted(glob(pathS))
    else:
        listFileS = [None for i in range(len(listFileU))]

    # fill temperature file list (list of None if computeH == False)
    if computeH:
        listFileT = sorted(glob(pathT))
    else:
        listFileT = [None for i in range(len(listFileU))]

    # check length of lists
    if len(listFileU) == 0:
        raise Exception("Error! list of U files is empty!")             
    if len(listFileU) != len(listFileV):
        raise Exception("Error! number of U and V files does not match! lenU: %d - lenV: %d" % (len(listFileU),len(listFileV)))     
    
    if computeS:
        if len(listFileS) != len(listFileU):
            raise Exception("Error! number of S and U files does not match! lenS: %d - lenU: %d" % (len(listFileS),len(listFileU)))

    if computeH:
        if len(listFileT) != len(listFileU):
            raise Exception("Error! number of T and U files does not match! lenT: %d - lenU: %d" % (len(listFileT),len(listFileU)))

    listFileZ = sorted(glob(pathZMesh))
    if len(listFileZ) == 1:
        print("vertical mesh is assumed constant")
        constZ = True

    # calculate static faces of cells along the spline for x and y points
    # yFaces = dy * dz ---> transp = yFaces * u (m**3/s)
    # xFaces = dz * dz ---> transp = xfaces * v (m**3/s)

    if transectType != "pureX":
        yFacesDa = calcYFaces(indicesU=indicesU,
               dsCoords=dsCoords,
               dsMeshZ=dsMeshZ,
               dsMaskU=dsMaskU,
               verboseLevel=verboseLevel)

    if transectType != "pureY":
        xFacesDa = calcXFaces(indicesV=indicesV,
               dsCoords=dsCoords,
               dsMeshZ=dsMeshZ,
               dsMaskV=dsMaskV,
               verboseLevel=verboseLevel)

    # calculate signs of segments depending on orientation
    # with respect to transect
    signsU, signsV = calc_signs(
                        indicesU = indicesU,
                        indicesV = indicesV,
                        latsF = latsF,
                        lonsF = lonsF,
                        latsTransect = latsTransect,
                        lonsTransect = lonsTransect,
                        dsMeshZ = dsMeshZ,
                        transectType = transectType,
                        verboseLevel = verboseLevel
                            )

    dates = []
    volumeTransportTotal = []
    volumeTransportDirection1 = []    
    volumeTransportDirection2 = []        
    volumeTransportZ = []
    saltTransportTotal = []
    heatTransportTotal = []
    
    if transectType == "general":
        for ifile,(fileU,fileV,fileS,fileT) in enumerate(zip(listFileU,listFileV,listFileS,listFileT)):  
            if verboseLevel>0:
                print("##########################")
                print("ifile: ",ifile)
                print("fileU: %s" % fileU)
                print("fileV: %s" % fileV)
            daU_sel = loadNEMOuVelocity(pathU=fileU,indicesU=indicesU)
            daV_sel = loadNEMOvVelocity(pathV=fileV,indicesV=indicesV)
            if verboseLevel>0:
                print(daU_sel.sizes)
                print(daV_sel.sizes)
            if verboseLevel>3:
                print(" ------ daU_sel ------")
                print(daU_sel)
                print(" ------ daV_sel ------")
                print(daV_sel)

            uTransp = yFacesDa.values * daU_sel.values * signsU # (Z,nUpoints)
            vTransp = xFacesDa.values * daV_sel.values * signsV # (Z,nVpoints)
            if verboseLevel>1:
                print("uTransp.shape:", uTransp.shape)
                print("vTransp.shape:", vTransp.shape)
            uTransp0    = np.nansum(uTransp) * volumeConversion        # 0 because now is 0-dimension
            vTransp0    = np.nansum(vTransp) * volumeConversion  
            volTransp0  = uTransp0 + vTransp0
            if verboseLevel>1:
                print("uTransp0: ", uTransp0)
                print("vTransp0: ", vTransp0)
            # direction 1 (>=0)
            uTransp0_direction1 = np.nansum(uTransp[uTransp>=0]) * volumeConversion
            vTransp0_direction1 = np.nansum(vTransp[vTransp>=0]) * volumeConversion
            volTransp0_dir1     = uTransp0_direction1 + vTransp0_direction1
            # direction 2 (<0)
            uTransp0_direction2 = np.nansum(uTransp[uTransp<0]) * volumeConversion    
            vTransp0_direction2 = np.nansum(vTransp[vTransp<0]) * volumeConversion
            volTransp0_dir2     = uTransp0_direction2 + vTransp0_direction2

            # horizontal sum of transport: sum_U (uTransp(z,nUpoints)) -->  uTransp(z)
            uTranspZ    = np.nansum(uTransp,axis=1) * volumeConversion
            vTranspZ    = np.nansum(vTransp,axis=1) * volumeConversion
            volTranspZ  = uTranspZ + vTranspZ        # (z)

            # store informations in lists
            dates.append(daU_sel["time_counter"].values)
            volumeTransportTotal.append(volTransp0)
            volumeTransportDirection1.append(volTransp0_dir1)
            volumeTransportDirection2.append(volTransp0_dir2)
            volumeTransportZ.append(volTranspZ)

            if computeS:
                print("  fileS: %s " % fileS)
                daTU_sel        = loadNEMOuTracer(pathT=fileS,varName=varNameS,indicesU=indicesU)
                daTV_sel        = loadNEMOvTracer(pathT=fileS,varName=varNameS,indicesV=indicesV)
                uSaltTransp     = uTransp * daTU_sel.values # (Z,nUpoints) 
                vSaltTransp     = vTransp * daTV_sel.values # (Z,nVpoints)
                uSaltTransp0    = rhoWater * np.nansum(uSaltTransp) * 1e-12     # to have units of 10^9 kg/s
                vSaltTransp0    = rhoWater * np.nansum(vSaltTransp) * 1e-12 
                saltTransp0     = uSaltTransp0 + vSaltTransp0
                saltTransportTotal.append(saltTransp0)

            if computeH:
                print("  fileT: %s " % fileS)
                daTU_sel        = loadNEMOuTracer(pathT=fileT,varName=varNameT,indicesU=indicesU)
                daTV_sel        = loadNEMOvTracer(pathT=fileT,varName=varNameT,indicesV=indicesV)
                uHeatTransp     = rhoWater * cpWater * uTransp * (daTU_sel.values - tempRef )   # (Z,nUpoints) 
                vHeatTransp     = rhoWater * cpWater * vTransp * (daTV_sel.values - tempRef )   # (Z,nVpoints)
                uHeatTransp0    = np.nansum(uHeatTransp) * 1e-15     # to have units of 10^15 Watts
                vHeatTransp0    = np.nansum(vHeatTransp) * 1e-15 
                heatTransp0     = uHeatTransp0 + vHeatTransp0
                heatTransportTotal.append(heatTransp0)

    if transectType == "pureX":
        for ifile,(fileU,fileV,fileS,fileT) in enumerate(zip(listFileU,listFileV,listFileS,listFileT)):  
            if verboseLevel>0:
                print("##########################")
                print("ifile: ",ifile)
                print("fileV: %s" % fileV)
            daV_sel = loadNEMOvVelocity(pathV=fileV,indicesV=indicesV)
            if verboseLevel>0:
                print(daV_sel.sizes)
            if verboseLevel>3:
                print(" ------ daV_sel ------")
                print(daV_sel)

            vTransp = xFacesDa.values * daV_sel.values * signsV # (Z,nVpoints)
            if verboseLevel>1:
                print("vTransp.shape:", vTransp.shape)
            vTransp0    = np.nansum(vTransp) * volumeConversion  
            volTransp0  =  vTransp0
            if verboseLevel>1:
                print("vTransp0: ", vTransp0)
            # direction 1 (>=0)
            vTransp0_direction1 = np.nansum(vTransp[vTransp>=0]) * volumeConversion
            volTransp0_dir1     = vTransp0_direction1
            # direction 2 (<0)
            vTransp0_direction2 = np.nansum(vTransp[vTransp<0]) * volumeConversion
            volTransp0_dir2     = vTransp0_direction2

            # horizontal sum of transport: sum_U (uTransp(z,nUpoints)) -->  uTransp(z)
            vTranspZ    = np.nansum(vTransp,axis=1) * volumeConversion
            volTranspZ  = vTranspZ        # (z)

            # store informations in lists
            dates.append(daU_sel["time_counter"].values)
            volumeTransportTotal.append(volTransp0)
            volumeTransportDirection1.append(volTransp0_dir1)
            volumeTransportDirection2.append(volTransp0_dir2)
            volumeTransportZ.append(volTranspZ)

            if computeS:
                print("  fileS: %s " % fileS)
                daTV_sel        = loadNEMOvTracer(pathT=fileS,varName=varNameS,indicesV=indicesV)
                vSaltTransp     = vTransp * daTV_sel.values # (Z,nVpoints)
                vSaltTransp0    = rhoWater * np.nansum(vSaltTransp) * 1e-12 
                saltTransp0     = vSaltTransp0
                saltTransportTotal.append(saltTransp0)

            if computeH:
                print("  fileT: %s " % fileS)
                daTV_sel        = loadNEMOvTracer(pathT=fileT,varName=varNameT,indicesV=indicesV)
                vHeatTransp     = rhoWater * cpWater * vTransp * (daTV_sel.values - tempRef )   # (Z,nVpoints)
                vHeatTransp0    = np.nansum(vHeatTransp) * 1e-15 
                heatTransp0     = vHeatTransp0
                heatTransportTotal.append(heatTransp0)

    if transectType == "pureY":
        # to be implemented
        for ifile,(fileU,fileV,fileS,fileT) in enumerate(zip(listFileU,listFileV,listFileS,listFileT)):  
            if verboseLevel>0:
                print("##########################")
                print("ifile: ",ifile)
                print("fileU: %s" % fileU)
            daU_sel = loadNEMOuVelocity(pathU=fileU,indicesU=indicesU)
            if verboseLevel>0:
                print(daU_sel.sizes)
            if verboseLevel>3:
                print(" ------ daU_sel ------")
                print(daU_sel)

            uTransp = yFacesDa.values * daU_sel.values * signsU # (Z,nUpoints)
            if verboseLevel>1:
                print("uTransp.shape:", uTransp.shape)
            uTransp0    = np.nansum(uTransp) * volumeConversion        # 0 because now is 0-dimension
            volTransp0  = uTransp0 
            if verboseLevel>1:
                print("uTransp0: ", uTransp0)
            # direction 1 (>=0)
            uTransp0_direction1 = np.nansum(uTransp[uTransp>=0]) * volumeConversion
            volTransp0_dir1     = uTransp0_direction1
            # direction 2 (<0)
            uTransp0_direction2 = np.nansum(uTransp[uTransp<0]) * volumeConversion    
            volTransp0_dir2     = uTransp0_direction2 

            # horizontal sum of transport: sum_U (uTransp(z,nUpoints)) -->  uTransp(z)
            uTranspZ    = np.nansum(uTransp,axis=1) * volumeConversion
            volTranspZ  = uTranspZ         # (z)

            # store informations in lists
            dates.append(daU_sel["time_counter"].values)
            volumeTransportTotal.append(volTransp0)
            volumeTransportDirection1.append(volTransp0_dir1)
            volumeTransportDirection2.append(volTransp0_dir2)
            volumeTransportZ.append(volTranspZ)

            if computeS:
                print("  fileS: %s " % fileS)
                daTU_sel        = loadNEMOuTracer(pathT=fileS,varName=varNameS,indicesU=indicesU)
                uSaltTransp     = uTransp * daTU_sel.values # (Z,nUpoints) 
                uSaltTransp0    = rhoWater * np.nansum(uSaltTransp) * 1e-12     # to have units of 10^9 kg/s
                saltTransp0     = uSaltTransp0
                saltTransportTotal.append(saltTransp0)

            if computeH:
                print("  fileT: %s " % fileS)
                daTU_sel        = loadNEMOuTracer(pathT=fileT,varName=varNameT,indicesU=indicesU)
                uHeatTransp     = rhoWater * cpWater * uTransp * (daTU_sel.values - tempRef )   # (Z,nUpoints) 
                uHeatTransp0    = np.nansum(uHeatTransp) * 1e-15     # to have units of 10^15 Watts
                heatTransp0     = uHeatTransp0 
                heatTransportTotal.append(heatTransp0)

    # create Output DataArray
    outFile     = outFileRoot % transectName
    datesArray  = np.array(dates) 
    depthArray  = yFacesDa["Z"].values

    # save Timeseries (time)
    daVolumeTransportTotal      = xr.DataArray(data = np.array(volumeTransportTotal),coords={"time":datesArray})
    daVolumeTransportTotal.name = "volume_transport_total"
    daVolumeTransportTotal      = daVolumeTransportTotal.assign_attrs({"units":"Sv","transect_name":transectName})
    
    daVolumeTransportDirection1 = xr.DataArray(data = np.array(volumeTransportDirection1),coords={"time":datesArray})
    daVolumeTransportDirection1.name = "volume_transport_direction1"
    daVolumeTransportDirection1 = daVolumeTransportDirection1.assign_attrs({"units":"Sv","transect_name":transectName})

    daVolumeTransportDirection2 = xr.DataArray(data = np.array(volumeTransportDirection2),coords={"time":datesArray})
    daVolumeTransportDirection2.name = "volume_transport_direction2"
    daVolumeTransportDirection2 = daVolumeTransportDirection2.assign_attrs({"units":"Sv","transect_name":transectName})

    # save horizontally integrated transport (time,Z)
    daVolumeTransportZ      = xr.DataArray( data = np.array(volumeTransportZ), coords = {"time":datesArray,"Z":depthArray})
    daVolumeTransportZ.name = "volume_transport_z"
    daVolumeTransportZ      = daVolumeTransportZ.assign_attrs({"units":"Sv","transect_name":transectName})

    if verboseLevel>0:
        print("Creating Output File")

    # merge arrays 
    # It's mandatory to assign names to DataArrays before merging!
    dsOut = xr.merge([daVolumeTransportTotal,
                      daVolumeTransportDirection1,
                      daVolumeTransportDirection2,
                      daVolumeTransportZ])

    if computeS:
        # create DataArray to store Salt transport
        daSaltTransportTotal = xr.DataArray( data = np.array(saltTransportTotal),coords={"time":datesArray} )
        daSaltTransportTotal.name = "salt_transport_total"
        daSaltTransportTotal = daSaltTransportTotal.assign_attrs({"units":"10^9 Kg/s","transect_name":transectName})
        # update dataset
        dsOut["salt_transport_total"] = daSaltTransportTotal

    if computeH:
        # create DataArray to store Salt transport
        daHeatTransportTotal = xr.DataArray( data = np.array(heatTransportTotal),coords={"time":datesArray} )
        daHeatTransportTotal.name = "heat_transport_total"
        daHeatTransportTotal = daHeatTransportTotal.assign_attrs({"units":"10^15 Watt","transect_name":transectName})
        # update dataset
        dsOut["heat_transport_total"] = daHeatTransportTotal

    # assign [m] units to Z coordinate
    dsOut.coords["Z"] = dsOut.coords["Z"].assign_attrs({"units":"m"})

    dsOut   = dsOut.assign_attrs({"transect_name":transectName})
    dsOut   = dsOut.assign_attrs({"transect_type":transectType})

    dsOut.to_netcdf(outFile,encoding={"time" : {"dtype": "i4"} },)

    #daVolumeTransportTotal.to_netcdf(outFile)

    #return indicesF,indicesU,indicesV


def calcXFaces(
                indicesV    = [],
                dsCoords    = xr.Dataset(),
                dsMeshZ     = xr.Dataset(),
                dsMaskV     = xr.Dataset(),
                verboseLevel    = 0
               ):

    # calculation of "X" faces (defined on V-points)

    #          Xface
    #     <-------------->    
    #
    #     .------v-------.     
    #     |              |   
    #     |              |   
    #     u      o       u       transport = v *  dx * dz 
    #     |              |                   v * e1v * e3v
    #     |              |   
    #     .------v-------.    
    #


    if verboseLevel >0:
        print("Calculating X Faces")

    # use extracted indices for U,V
    tgt_y_v     = xr.DataArray(np.array([idx[0] for idx in indicesV]),dims="vpoints")
    tgt_x_v     = xr.DataArray(np.array([idx[1] for idx in indicesV]),dims="vpoints")
    e1v_sel     = dsCoords["e1v"].isel(x=tgt_x_v,y=tgt_y_v).squeeze()
    e3v_sel     = dsMeshZ["e3v_0"].isel(X=tgt_x_v,Y=tgt_y_v).squeeze()              #(Z,nUpoints)
    maskV_sel   = dsMaskV["mask"].isel(x=tgt_x_v,y=tgt_y_v).squeeze()               #(depthu,nVpoints)

    # expand dims of e2u (nu) -> (z,nu)
    nz = e3v_sel.sizes["Z"]
    e1v_sel = e1v_sel.expand_dims(dim={"Z":nz})

    maskV_sel = maskV_sel.rename({"depthv":"Z"})

    if verboseLevel > 0:
        print("-- e1v_sel --")        
        print(e1v_sel.sizes)
        print("-- e3v_sel --")
        print(e3v_sel.sizes)                
        print("-- maskV_sel --")
        print(maskV_sel.sizes)

    # use numpy arrays for the multiplication
    xFaces = e1v_sel.values * e3v_sel.values * maskV_sel.values

    if len(indicesV) == 1:                          # (z,) -> (z,nVpoints)
        xFaces = np.expand_dims(xFaces,axis=1)

    xFacesDa = xr.DataArray(xFaces,dims=["Z","vpoints"])

    # add informations on depth
    xFacesDa = xFacesDa.assign_coords({"Z":maskV_sel["Z"]})

    if verboseLevel > 1:
        print("---- xFaces shape ----")        
        print(xFacesDa.sizes)

    return xFacesDa

def calcYFaces(
                indicesU    = [],
                dsCoords    = xr.Dataset(),
                dsMeshZ     = xr.Dataset(),
                dsMaskU     = xr.Dataset(),
                verboseLevel     = False
               ):

    # calculation of "Y" faces (defined on U-points)

    #    
    #     .------v-------.   ^  
    #     |              |   |
    #     |              |   |
    #     u      o       u   |  Yface      transport = u *  dy * dz 
    #     |              |   |                         u * e2u * e3u
    #     |              |   |
    #     .------v-------.   v 
    #

    if verboseLevel>0:
        print("Calculating Y Faces")

    tgt_y_u     = xr.DataArray(np.array([idx[0] for idx in indicesU]),dims="upoints")
    tgt_x_u     = xr.DataArray(np.array([idx[1] for idx in indicesU]),dims="upoints")
    # extract Arrays 
    e2u_sel     = dsCoords["e2u"].isel(x=tgt_x_u,y=tgt_y_u).squeeze()               #(nUpoints)
    e3u_sel     = dsMeshZ["e3u_0"].isel(X=tgt_x_u,Y=tgt_y_u).squeeze()              #(Z,nUpoints)
    maskU_sel   = dsMaskU["mask"].isel(x=tgt_x_u,y=tgt_y_u).squeeze()               #(depthu,nUpoints)

    # expand dims of e2u (nu) -> (z,nu)
    nz = e3u_sel.sizes["Z"]
    e2u_sel = e2u_sel.expand_dims(dim={"Z":nz})

    maskU_sel = maskU_sel.rename({"depthu":"Z"})

    if verboseLevel>1:
        print("-- e2u_sel --")        
        print(e2u_sel.sizes)
        print("-- e3u_sel --")
        print(e3u_sel.sizes)                
        print("-- maskU_sel --")
        print(maskU_sel.sizes)

    # use numpy arrays for the multiplication
    yFaces = e2u_sel.values * e3u_sel.values * maskU_sel.values

    if len(indicesU) == 1:                          #(z,) -> (z,nUpoints)
        yFaces = np.expand_dims(yFaces,axis=1)

    yFacesDa = xr.DataArray(yFaces,dims=["Z","upoints"])

    # add informations on depth
    yFacesDa = yFacesDa.assign_coords({"Z":maskU_sel["Z"]})

    if verboseLevel>1:
        print("---- yFaces shape ----")        
        print(yFacesDa.sizes)

    return yFacesDa


def calcNEMOTransportSingle(                            # calculate the transport for a single Time-Level
                        indicesV = [],
                        indicesU = [],                              
                        dsCoords = xr.Dataset(),        # xr.Dataset of 2D coordinates.nc
                        dsZMesh = xr.Dataset(),         # xr.Dataset of file with vertical scale factors   
                        daU = xr.DataArray(),           # xr.DataArray with U velocity
                        daV = xr.DataArray()            # xr.DataArray with V velocity                                   
                    ):
    
    # stote V indices in DataArrays
    tgt_y_v     = xr.DataArray(np.array([idx[0] for idx in indicesV]),dims="vpoints")
    tgt_x_v     = xr.DataArray(np.array([idx[1] for idx in indicesV]),dims="vpoints")

    # stote U indices in DataArrays
    tgt_y_u     = xr.DataArray(np.array([idx[0] for idx in indicesU]),dims="upoints")
    tgt_x_u     = xr.DataArray(np.array([idx[1] for idx in indicesU]),dims="upoints")



    pass

def loadNEMOuVelocity(
                    pathU = "",
                    indicesU = [],
                    varName = "uo"                    
                    ):

    # store indices in DataArrays      
    tgt_y_u     = xr.DataArray(np.array([idx[0] for idx in indicesU]),dims="upoints")
    tgt_x_u     = xr.DataArray(np.array([idx[1] for idx in indicesU]),dims="upoints")

    if not isfile(pathU):
        raise Exception("Error! file %s not existing" % pathU)
    
    dsU = xr.open_dataset(pathU)
    daU = dsU[varName]
    timeStamp = dsU["time_counter"]
    # make sure I don't loose time info
    daU = daU.assign_coords({"time_counter":timeStamp})

    daU_sel     = daU.isel(x=tgt_x_u  ,y=tgt_y_u).squeeze()     # (Z,nUpoints)       # note that NEMO has x,y not X,Y

    # particular case of 1 single value u point value
    if len(indicesU) == 1:
        daU_sel = daU_sel.expand_dims(dim={"upoints":1},axis=1)

    dsU.close()

    return daU_sel

def loadNEMOvVelocity(
                    pathV = "",
                    indicesV = [],
                    varName = "vo"
                    ):

    # access to a single dataset and get the DataArray of U-velocity

    # store indices in DataArrays      
    tgt_y_v     = xr.DataArray(np.array([idx[0] for idx in indicesV]),dims="vpoints")
    tgt_x_v     = xr.DataArray(np.array([idx[1] for idx in indicesV]),dims="vpoints")

    if not isfile(pathV):
        raise Exception("Error! file %s not existing" % pathV)
    
    dsV = xr.open_dataset(pathV)
    daV = dsV[varName]
    timeStamp = dsV["time_counter"]
    # make sure I don't loose time info
    daV = daV.assign_coords({"time_counter":timeStamp})

    daV_sel = daV.isel(x=tgt_x_v,y=tgt_y_v).squeeze()   # (Z,nVpoints)      # note that NEMO has x,y not X,Y 

    if len(indicesV) == 1:
        daV_sel = daV_sel.expand_dims(dim={"vpoints":1},axis=1)

    dsV.close()

    return daV_sel

def loadNEMOuTracer(
                    pathT = "",
                    varName = "",
                    indicesU = [],
                    daU_sel = xr.DataArray()
                    ):

    #       Estimate tracer on U points with average
    #       
    #       If one of i or i+1 tracer points is 0 the average is not accurate
    #       but there is no effect on the transport calculation, since
    #       also the u is 0.      
    #
    #      u      T      u      T      u      T      u
    #      |      o      |      o      |      o      |
    #            i-1            i            i+1                       
    #     i-2           i-1            i            i+1
    #

    # store indices in DataArrays      
    tgt_y_u     = xr.DataArray(np.array([idx[0] for idx in indicesU]),dims="upoints")
    tgt_x_u     = xr.DataArray(np.array([idx[1] for idx in indicesU]),dims="upoints")

    if not isfile(pathT):
        raise Exception("Error! file %s not existing" % pathT)

    dsT = xr.open_dataset(pathT)
    daT = dsT[varName]
    timeStamp = dsT["time_counter"]
    # make sure I don't loose time info
    daT = daT.assign_coords({"time_counter":timeStamp})

    daTU_sel_ip1 = daT.isel(x=tgt_x_u+1,y=tgt_y_u).squeeze()     # (Z,nUpoints) 
    daTU_sel     = daT.isel(x=tgt_x_u  ,y=tgt_y_u).squeeze()     # (Z,nUpoints) 

    """
    # loop on indicesU to select upwind tracer point

    nx = daT.sizes["X"]
    # loop on z
    for z in range(outArray.shape[0]):
        # loop on horizontal indices
        for ipoint,(j,i) in enumerate(indicesU):
            uvel = daU_sel.squeeze().isel(Z=z,upoints=ipoint)
            # select upwind index
            xidx = int(i+max(0,-np.sign(uvel.item())))
            # make sure I don't go out of boundaries
            xidx = np.clip(xidx,0,nx-1)
            outArray[z,ipoint] = daT.isel(Z=z,Y=j,X=xidx).squeeze().item()

    """

    # replace data
    daTU_sel.data = 0.5 * (daTU_sel.data + daTU_sel_ip1.data)

    dsT.close()

    return daTU_sel

def loadNEMOvTracer(
                    pathT = "",
                    varName = "",
                    indicesV = [],
                    daV_sel = xr.DataArray()
                    ):

    #       Estimate tracer on V points with average
    #       
    #       If one of j or j+1 tracer points is 0 the average is not accurate
    #       but there is no effect on the transport calculation, since
    #       also the V is 0.

    #
    #      v  ---     j+1
    # 
    #      T   o  j+1
    #
    #      v  ---     j
    #
    #      T   o  j   
    #
    #      v  ---     j-1
    #  
    #      T   o  j-1
    #
    #      v  ---     j-2
    #

    # store indices in DataArrays      
    tgt_y_v     = xr.DataArray(np.array([idx[0] for idx in indicesV]),dims="vpoints")
    tgt_x_v     = xr.DataArray(np.array([idx[1] for idx in indicesV]),dims="vpoints")

    if not isfile(pathT):
        raise Exception("Error! file %s not existing" % pathT)

    dsT = xr.open_dataset(pathT)
    daT = dsT[varName]
    timeStamp = dsT["time_counter"]
    # make sure I don't loose time info
    daT = daT.assign_coords({"time_counter":timeStamp})

    daTV_sel_jp1 = daT.isel(x=tgt_x_v  ,y=tgt_y_v+1).squeeze()     # (Z,nVpoints) 
    daTV_sel     = daT.isel(x=tgt_x_v  ,y=tgt_y_v  ).squeeze()     # (Z,nVpoints) 

    """
    # loop on indicesV to select upwind tracer point
    # upwind search is really expensive

    outArray = np.zeros((daV_sel.shape))      # shape of outArray is (Z,)     
    
    ny = daT.sizes["Y"]
    # loop on Z
    for z in range(outArray.shape[0]):
        # loop on horizontal indices
        for ipoint,(j,i) in enumerate(indicesV):
            vvel = daV_sel.squeeze().isel(Z=z,vpoints=ipoint)
            # select upwind index
            yidx = int(j+max(0,-np.sign(vvel.item())))
            # make sure I don't go out of boundaries
            yidx = np.clip(yidx,0,ny-1)
            outArray[z,ipoint] = daT.isel(Z=z,Y=yidx,X=i).squeeze().item()
    """
            
    # replace data
    daTV_sel.data = 0.5 * (daTV_sel.data + daTV_sel_jp1.data)

    dsT.close()

    return daTV_sel
