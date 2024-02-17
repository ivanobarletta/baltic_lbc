# coding: utf-8
varSiren = dsSiren[varName].isel(T=0,Z=iLevel)
varCDO = dsCDO[varName].isel(T=0,Z=iLevel) 
fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(8,5))
varSiren.plot(ax=ax[0],cmap="jet",x="nav_lon",y="nav_lat",vmin=0,vmax=6)
varCDO.plot(ax=ax[1],cmap="jet",x="nav_lon",y="nav_lat",vmin=0,vmax=6)
