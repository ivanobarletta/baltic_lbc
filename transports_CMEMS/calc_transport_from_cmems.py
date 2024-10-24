import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np

"""
    Script to calculate zonal transport (meridional has yet to be implemented) 
    from a CMEMS product

    1) Necessary files

        Make sure you have the static files containing

        1a) horizontal scale factors (set pathE2T)
        1b) verttical scale factors (set pathE3T)

        the files containing the variables needed for the 
        calculation of transports

        1c) the file with u-velocity

        IF you set computeS to True (calculation of salinity transport)

            1d) the path of file containing salinity (pathS) 
                and varnameS

        IF you set computeH to True (calculation of heat transport)

            1e) the path of file containing temperature (pathT) 
                and varnameT
        

    2) NOTE WELL!

        the script assumes that the spatial domain
        of all files (static, variables) is the same!

    3) Transect info

        for vertical transect you must set

        transectCoords = [lon,lat1,lat2]


"""



def calc_x_transport_single_file(
                            parentModel="model",
                            transectName="name",
                            transectCoords = [10,50,60],                            
                            pathU="",
                            pathS=None,
                            pathT=None,
                            pathSSH="",
                            pathE2T="",
                            pathE3T="",
                            pathH="",
                            uVarName = "uo",
                            sshVarName = "zos_detided",
                            saltVarName = "so",
                            thetaVarName = "thetao",
                            depthVarName = "deptho",
                            dates = None
                            ):
    
    # the variables needed are all in a single file for each variable (time,depth,lat,lon)

    tref = 0 

    cp          = 3996              # J / kg / K    specific heat to sea water
    rho         = 1026              # kg / m**3     water density
    tref        = 0

    volumeConversion    = 1e-6      # multiply to get units of 10^6 m**3/s (Sv)    
    saltConversion      = 1e-12     # multiply to get units of 10^9 Kg/s    
                                    # rho * s * volumeTransport
                                    # (kg / m^3) * (g / kg) * (m^3 / s) = g/s
                                    # I need a factor of 10^(-12) to get units of 10^9 Kg/s    
    heatConversion      = 1e-15     # multiply to get units of 10^15 W

    xTransect = transectCoords[0]
    yTransect1 = transectCoords[1]
    yTransect2 = transectCoords[2]

    if pathS != None:
        so      = xr.open_dataset(pathS)[saltVarName]       # (t,z,y,x)
    if pathT != None:
        thetao  = xr.open_dataset(pathT)[thetaVarName]      # (t,z,y,x)        
    uo  = xr.open_dataset(pathU)[uVarName]                  # (t,z,y,x)
    ssh = xr.open_dataset(pathSSH)[sshVarName]              # (t,y,x)
    e2t = xr.open_dataset(pathE2T)["e2t"]                   # (y,x)
    e3t = xr.open_dataset(pathE3T)["e3t"]                   # (z,y,x)
    H   = xr.open_dataset(pathH)[depthVarName]              # (y,x)

    # reduce all dataarrays
    print("slicing dataarrays")
    if pathS != None:
        so  = so.sel(longitude=xTransect,method="nearest").sel(latitude=slice(yTransect1,yTransect2))           #(t,z,y)
        if dates != None:
            print("time slice")
            so = so.sel(time=slice(dates[0],dates[1]))
    if pathT != None:
        if dates != None:
            print("time slice")
            thetao = thetao.sel(time=slice(dates[0],dates[1]))
        thetao  = thetao.sel(longitude=xTransect,method="nearest").sel(latitude=slice(yTransect1,yTransect2))   #(t,z,y)

    uo  = uo.sel(longitude=[xTransect],method="nearest").sel(latitude=slice(yTransect1,yTransect2))   #(t,z,y)
    ssh = ssh.sel(longitude=[xTransect],method="nearest").sel(latitude=slice(yTransect1,yTransect2))  #(t,y)
    e2t = e2t.sel(longitude=[xTransect],method="nearest").sel(latitude=slice(yTransect1,yTransect2))  #(y)
    e3t = e3t.sel(longitude=[xTransect],method="nearest").sel(latitude=slice(yTransect1,yTransect2))  #(z,y)
    H   = H.sel(longitude=[xTransect],method="nearest").sel(latitude=slice(yTransect1,yTransect2))    #(y)

    print(" uo.sizes ", uo.sizes)
    print("ssh.sizes ", ssh.sizes)
    print("e2t.sizes ", e2t.sizes)
    print("e3t.sizes ", e3t.sizes)
    print("  H.sizes ", H.sizes)

    # check on depth size of uo and e3t (possibly different)
    if uo.sizes["depth"] != e3t.sizes["depth"]:
        print("uo and e3t have different nz. Fixing..")
        nz_max = max(uo.sizes["depth"],e3t.sizes["depth"])
        uo = uo.sel(depth=slice(None,nz_max))
        e3t = e3t.sel(depth=slice(None,nz_max))
        if pathS != None:
            so = so.sel(depth=slice(None,nz_max))
        if pathT != None:
            thetao = thetao.sel(depth=slice(None,nz_max))

    if dates != None:
        uo  = uo.sel(time=slice(dates[0],dates[1]))
        ssh = ssh.sel(time=slice(dates[0],dates[1]))        

    # H(lat,lon) -> H(time,lat,lon)
    H = H.expand_dims(dim={"time":ssh.coords["time"]})  #(t,y)
    
    # calculate Jacobian (xarray complains if coordinates are not exactly the same!)
    H["longitude"] = ssh["longitude"]
    H["latitude"] = ssh["latitude"]

    Jacobian = (H+ssh) / H      #(t,y)

    Jacobian = Jacobian.expand_dims(dim={"depth":e3t.coords["depth"]},axis=1)   # (t,z,y)

    e2t["latitude"] = Jacobian["latitude"]
    e3t["latitude"] = Jacobian["latitude"]
    e2t["longitude"] = Jacobian["longitude"]
    e3t["longitude"] = Jacobian["longitude"]

    # This works also If I don't expand dimension of e2t/e3t
    # N.B. faces are not masked. I rely on the masked values of either U or S to 
    # compute transport properly
    # Note. (2024-10-15) I had to add the brackets. It seems that the order of operations
    # is important
    
    xfaces = Jacobian * ( e3t * e2t )                                               # (time,depth,latitude,longitude) (longitude=1..)

    # masking of xfaces ( I do this because scale factors are not masked for BALMFC while 
    # they are for GLOMFC)

    mask3D = xr.full_like(uo.isel(time=0),fill_value=np.nan)
    mask3D.data[~np.isnan(uo.isel(time=0).data)] = 1
    xfaces2  = xfaces * mask3D 

    print("xfaces shape", xfaces.shape)

    uo["latitude"]  = xfaces["latitude"]
    uo["depth"]     = xfaces["depth"]

    print("uo shape", uo.shape)

    print("calculating Volume Transports")
    volumeTransport_data    = uo.data * xfaces.data                                    # (time,depth,latitude,longitude) [m^3/s]
    volumeTransport         = xr.full_like(uo,fill_value=np.nan)
    volumeTransport.data    = volumeTransport_data * volumeConversion
    volumeTransport.name    = "volume_transport2D"
    volumeTransport         = volumeTransport.assign_attrs({"units":"Sv"})

    # separate direction1 from direction2
    print("calculating Volume Transports_direction1")
    volumeTransport_direction1 = volumeTransport.where(volumeTransport>=0)  # (time,depth,latitude) [m^3/s]
    volumeTransport_direction2 = volumeTransport.where(volumeTransport<0)   # (time,depth,latitude) [m^3/s]

    outFile = "transports_%s_%s.nc" % (parentModel, transectName)

    dsOut   = xr.Dataset()

    xfaces2.name = "xfaces"
    xfaces2  = xfaces2.assign_attrs({"units":"m2"})
    # this appears to be important not to spoil the outcome
    # If I don't do this, the transport data appear all NaN for the GLOMFC (I don't know why...)
    xfaces2["time"] = volumeTransport["time"]       

    dsOut["xfaces"] = xfaces2

    dsOut["volume_transport2D"] = volumeTransport



    print("filling dataset with timeseries")
    # fill dataset with timeseries 
    volumeTransportTimeseries = volumeTransport.sum(dim=["depth","latitude"])
    #volumeTransportTimeseries *= volumeConversion
    volumeTransportTimeseries = volumeTransportTimeseries.assign_attrs({"units":"Sv"})
    dsOut["volumeTransport_total"] = volumeTransportTimeseries

    print("filling dataset with timeseries")
    volumeTransportTimeseries = volumeTransport_direction1.sum(dim=["depth","latitude"])
    #volumeTransportTimeseries *= volumeConversion
    volumeTransportTimeseries = volumeTransportTimeseries.assign_attrs({"units":"Sv"})
    dsOut["volumeTransport_direction1"] = volumeTransportTimeseries

    print("filling dataset with timeseries")
    volumeTransportTimeseries = volumeTransport_direction2.sum(dim=["depth","latitude"])
    #volumeTransportTimeseries *= volumeConversion
    volumeTransportTimeseries = volumeTransportTimeseries.assign_attrs({"units":"Sv"})
    dsOut["volumeTransport_direction2"] = volumeTransportTimeseries

    volumeTransportHovmoller    = volumeTransport.sum(dim="latitude")
    #volumeTransportHovmoller    *= volumeConversion
    volumeTransportHovmoller    = volumeTransportHovmoller.assign_attrs({"units":"Sv"})
    dsOut["volumeTransportZ"]   = volumeTransportHovmoller


    if pathS != None:
        print("calculating Salt Transports")    
        so["latitude"]  = xfaces["latitude"]
        so["depth"]     = xfaces["depth"]
        saltTransport   = rho * so * volumeTransport * (1 / volumeConversion)                   # [g/s]
        saltTransportTimeseries = saltTransport.sum(dim=["depth","latitude"])
        saltTransportTimeseries *= saltConversion
        saltTransportTimeseries = saltTransportTimeseries.assign_attrs({"units":"10^9 Kg/s"})
        dsOut["saltTransport_total"] = saltTransportTimeseries
        saltTransport_direction1   = rho * so * volumeTransport_direction1                 # [g/s]
        saltTransportTimeseries = saltTransport_direction1.sum(dim=["depth","latitude"])
        saltTransportTimeseries *= saltConversion                                           
        saltTransportTimeseries = saltTransportTimeseries.assign_attrs({"units":"10^9 Kg/s"})
        dsOut["saltTransport_direction1"] = saltTransportTimeseries
        saltTransport_direction2   = rho * so * volumeTransport_direction2                 # [g/s]
        saltTransportTimeseries = saltTransport_direction2.sum(dim=["depth","latitude"])
        saltTransportTimeseries *= saltConversion                                           
        saltTransportTimeseries = saltTransportTimeseries.assign_attrs({"units":"10^9 Kg/s"})
        dsOut["saltTransport_direction2"] = saltTransportTimeseries

    if pathT != None:
        print("calculating Heat Transports")    
        thetao["latitude"]  = xfaces["latitude"]
        thetao["depth"]     = xfaces["depth"]
        heatTransport   = cp * (thetao - tref) * volumeTransport * (1/volumeConversion)       # [J/s] (Watts)
        heatTransportTimeseries = heatTransport.sum(dim=["depth","latitude"])
        heatTransportTimeseries *= heatConversion
        heatTransportTimeseries = heatTransportTimeseries.assign_attrs({"units":"10^15 W"})
        dsOut["heatTransport_total"] = heatTransportTimeseries
        heatTransport_direction1   = cp * (thetao - tref) * volumeTransport_direction1      # [J/s] (Watts)
        heatTransportTimeseries = heatTransport_direction1.sum(dim=["depth","latitude"])
        heatTransportTimeseries *= heatConversion
        heatTransportTimeseries = heatTransportTimeseries.assign_attrs({"units":"10^15 W"})
        dsOut["heatTransport_direction1"] = heatTransportTimeseries
        heatTransport_direction2   = cp * (thetao - tref) * volumeTransport_direction2      # [J/s] (Watts)
        heatTransportTimeseries = heatTransport_direction2.sum(dim=["depth","latitude"])
        heatTransportTimeseries *= heatConversion
        heatTransportTimeseries = heatTransportTimeseries.assign_attrs({"units":"10^15 W"})
        dsOut["heatTransport_direction2"] = heatTransportTimeseries

    dsOut   = dsOut.assign_attrs({"transectName":transectName,"transectLon":transectCoords[0],"transectLat1":transectCoords[1],"transectLat2":transectCoords[2]})

    print("saving to netcdf")
    dsOut.to_netcdf(outFile)


# note that the same transect must have slightly 
# different coordinates for the different parent models
# They have very different resolutions and mask
# and I might not have much sense to provide exactly the same transect 
# coordinates to both parent models

"""
parentModel     = "balmfc"

transectCoords  = [13.60,54.54,55.45]
transectName    = "arkona2"

transectCoords  = [13.37,54.65,55.45]
transectName    = "arkona"

root = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/BALMFC_PROD/"

pathU   = root + "cmems_mod_bal_phy_anfc_P1D-m_uo_9.04E-13.99E_53.01N-59.99N_0.50-591.94m_2022-01-01-2024-01-01.nc"             # u-velocity
pathS   = root + "cmems_mod_bal_phy_anfc_P1D-m_so_9.04E-13.99E_53.01N-59.99N_0.50-591.94m_2022-01-01-2024-01-01.nc"             # salinity
pathSSH = root + "cmems_mod_bal_phy-ssh_anfc_detided_P1D-m_zos_detided_9.04E-13.99E_53.01N-59.99N_2022-01-01-2024-01-01.nc"     # ssh
pathE2T = root + "cmems_mod_bal_phy_anfc_static_e1t-e2t-e3t_9.04E-13.99E_53.01N-59.99N_0.50-591.94m.nc"                         # dy
pathE3T = root + "cmems_mod_bal_phy_anfc_static_e1t-e2t-e3t_9.04E-13.99E_53.01N-59.99N_0.50-591.94m.nc"                         # dz
pathH   = root + "cmems_mod_bal_phy_anfc_static_multi-vars_9.04E-13.99E_53.01N-59.99N_0.50-591.94m.nc"                          # bathymetry
pathT   = root + "cmems_mod_bal_phy_anfc_P1D-m_thetao_9.04E-13.99E_53.01N-59.99N_0.50-591.94m_2022-01-01-2024-01-01.nc"         # temperature
sshVarName = "zos_detided"

"""

parentModel     = "glomfc"



transectCoords  = [13.58,54.5,55.5]
transectName    = "arkona2"

transectCoords  = [13.33,54.35,55.5]
transectName    = "arkona"

root = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/GLOMFC_PROD/"

pathU   = root + "DAILY/glo_anfc_subset.nc"             # u-velocity
pathS   = root + "DAILY/glo_anfc_subset.nc"             # salinity
pathSSH = root + "SSH/zos_subset.nc"                    # ssh
pathE2T = root + "static_box_12_15_53_56/cmems_mod_glo_phy_anfc_0.083deg_static_e1t-e2t-e3t_9.00E-14.00E_52.00N-60.00N_0.49-541.09m.nc"                         # dy
pathE3T = root + "static_box_12_15_53_56/cmems_mod_glo_phy_anfc_0.083deg_static_e1t-e2t-e3t_9.00E-14.00E_52.00N-60.00N_0.49-541.09m.nc"                         # dz
pathH   = root + "static_box_12_15_53_56/cmems_mod_glo_phy_anfc_0.083deg_static_multi-vars_9.00E-14.00E_52.00N-60.00N_0.49-541.09m.nc"                          # bathymetry
pathT   = root + "DAILY/glo_anfc_subset.nc"             # temperature
sshVarName = "zos"




dates   = ("2021-12-29","2023-12-31")
#dates   = ("2021-12-29","2022-4-30")

calc_x_transport_single_file(
                            parentModel=parentModel,
                            transectName=transectName,
                            transectCoords = transectCoords,                            
                            pathU=pathU,
                            pathS=pathS,
                            pathT=pathT,
                            pathSSH=pathSSH,
                            pathE2T=pathE2T,
                            pathE3T=pathE3T,
                            pathH=pathH,
                            dates=dates,
                            sshVarName=sshVarName
                            )

