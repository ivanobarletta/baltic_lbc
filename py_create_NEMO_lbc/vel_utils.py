from os.path import isfile
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

def project_velocity(ds,nml):

    # if ln_proj_vel is True project
    # geographical velocity to NEMO (i,j)
    # curvilinear system
    lnProJVel = nml["namvar"]["ln_proj_vel"]

    if lnProJVel:
        print("you want to project the velocity! Good Luck!")

        # get coordinates of F-points     
        lonsf, latsf = get_target_f_coords(nml)

        print("dims of lonsf,latsf (Must be +1 than T-points in each dimension)")
        print(lonsf.shape)
        print(latsf.shape)

        # build normals to cell edges
        iNorm_x,iNorm_y,jNorm_x,jNorm_y = build_ij_normal(lonsF=lonsf,latsF=latsf)
        
        # get vozocrtx,vomecrty
        vozo = ds.vozocrtx.T
        vome = ds.vomecrty.T
        
        print("iNorm_x.shape",iNorm_x.shape)
        print("iNorm_y.shape",iNorm_y.shape)
        print("jNorm_x.shape",jNorm_x.shape)
        print("jNorm_y.shape",jNorm_y.shape)        
        print("vozo.shape",vozo.shape)
        print("vome.shape",vome.shape)

        #scal_i = vozo.isel(T=0,Z=0).values * iNorm_x + vome.isel(T=0,Z=0).values * iNorm_y
        #print("scal_i.shape",scal_i.shape)

        print("Projecting the Velocity Field...")
        # do the projection       
        for t in range(ds.dims["T"]):
            print("    T = ", t)
            for z in range(ds.dims["Z"]):
                print ("      Z = ", z)     
                # make scalar product
                scal_i = vozo.isel(T=t,Z=z).values * iNorm_x + vome.isel(T=t,Z=z).values * iNorm_y
                scal_j = vozo.isel(T=t,Z=z).values * jNorm_x + vome.isel(T=t,Z=z).values * jNorm_y
                vozo.isel(T=t,Z=z).values[:] = scal_i[:]
                vome.isel(T=t,Z=z).values[:] = scal_j[:]        


        # update dataset
        ds["vozocrtx"] = vozo
        ds["vomecrty"] = vome

    return ds

#-----------------------------------------------------

def get_target_f_coords(nml):
    # get coordinates of F-points 
    #       
    #    F    V    F    V     F
    #
    #    U    T    U    T     U      
    #
    #    F    V    F    V     F
    #
    #    U    T    U    T     U      
    #
    #    F    V    F    V     F
    #  
    #  NxM T grid --> (N+1)x(M+1) F grid

    path = nml["namtgt"]["cn_coords_domain"]
    print("Searching for: %s" % path)
    if not isfile(path):
        print("the path of cn_coords_domain not present")
        print("the program needs F-points coordinates")
        raise Exception("Problem with cn_coords_domain")

    try:
        ds = xr.open_dataset(path)
    except:
        ds = xr.open_dataset(path,decode_times=False)    

    ds = ds.isel(time=0,z=0)

    lonsF = ds.glamf
    latsF = ds.gphif

    lonsT = ds.nav_lon
    latsT = ds.nav_lat
    ds.close()

    idxList = nml["nambdy"]["in_idx_list"]

    print ("idxList:")
    print (idxList)

    #               1 to 0 based numbering    
    #                  v 
    xIdx1 = idxList[0]-1
    xIdx2 = idxList[1]-1
    yIdx1 = idxList[2]-1
    yIdx2 = idxList[3]-1

    # subset coordinates F
    lonsF = lonsF.isel(x=slice(xIdx1-1,xIdx2+1),y=slice(yIdx1-1,yIdx2+1))
    latsF = latsF.isel(x=slice(xIdx1-1,xIdx2+1),y=slice(yIdx1-1,yIdx2+1))

    # subset coordinates T 
    lonsT = lonsT.isel(x=slice(xIdx1  ,xIdx2+1),y=slice(yIdx1  ,yIdx2+1))
    latsT = latsT.isel(x=slice(xIdx1  ,xIdx2+1),y=slice(yIdx1  ,yIdx2+1))

    """    
    plt.figure()
    plt.scatter(ds.nav_lon,ds.nav_lat,marker="o",color="k",s=0.4)
    plt.scatter(lonsT,latsT,marker="x",color="r",s=0.4)
    plt.scatter(lonsF,latsF,marker="v",color="g",s=0.2)
    plt.plot(lonsF,latsF,color="g")
    plt.plot(lonsF.T,latsF.T,color="g")
    plt.xlim(13,14.5)
    plt.ylim(54,56)
    plt.show()
    """

    return lonsF,latsF

#-----------------------------------------------------

def build_ij_normal(lonsF,latsF):
    # build normal vector to F grid 

    #                          j 
    #               F____     /       
    #              /      ------ ______F
    #             /                   /      
    #            /         T         /       
    #           /                   /  -> i       
    #          F ____              /
    #                 ----- ______F
    #    
    #    

    # vector pointing from
    # (i+1,j) to (i+1,j+1) called DJ    
    vDJ_x = lonsF.sel({"y": slice(1,None),"x":slice(1,None)}) - lonsF.sel({"y": slice(0,-1),"x":slice(1,None)})
    vDJ_y = latsF.sel({"y": slice(1,None),"x":slice(1,None)}) - latsF.sel({"y": slice(0,-1),"x":slice(1,None)})

    # vector pointing from 
    # (i,j+1) to (i+1,j+1) called DI        
    vDI_x = lonsF.sel({"x": slice(1,None),"y":slice(1,None)}) - lonsF.sel({"x": slice(0,-1),"y":slice(1,None)})
    vDI_y = latsF.sel({"x": slice(1,None),"y":slice(1,None)}) - latsF.sel({"x": slice(0,-1),"y":slice(1,None)})

    # rotate (vDJ_x,vDJ_y) 90 deg clockwise
    iNorm_x =   vDJ_y
    iNorm_y = - vDJ_x

    # rotate (vDI_x,vDI_y) 90 deg anti-clockwise
    jNorm_x = - vDI_y
    jNorm_y =   vDI_x
    
    # do some checks
    # dot ( vDI * vDJ  ) should be 0
    # dot ( vDJ * iNorm) should be 0
    # dot ( vDI * jNorm) should be 0

    """
    dot0 = vDI_x * vDJ_x + vDI_y * vDJ_y

    dot1 = vDJ_x * iNorm_x + vDJ_y * iNorm_y

    dot2 = vDI_x * jNorm_x + vDI_y * jNorm_y 

    print("testing the normal vectors")
    print()

    print(dot0)
    print()
    print(dot1)
    print()
    print(dot2)
    """

    


    # Normalize (iNorm,jNorm)
    iNormMod = np.sqrt(iNorm_x * iNorm_x + iNorm_y * iNorm_y)
    jNormMod = np.sqrt(jNorm_x * jNorm_x + jNorm_y * jNorm_y)

    iNorm_x /= iNormMod
    iNorm_y /= iNormMod
    jNorm_x /= jNormMod
    jNorm_y /= jNormMod

    return iNorm_x,iNorm_y,jNorm_x,jNorm_y     
