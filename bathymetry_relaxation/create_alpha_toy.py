import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 


# I want to create a function like this
#
# -------------------------__
#                            --__       
#                                ---_______   
# o   o   o   o   o   o   o   o   o   o   o       
#
# 0   1   2   3   4   5   6   7   8   9   10  
#
# f(6) = 1
# f(9) = 0
# m = (0 - 1) / (9 - 6) = -1/3
#
# f(i) = -1/3 * i + q
# 
# In general 
# m = -1 / (idx2-idx1)
# 
# f(i)

nx = 1093


idx1 = 1078-1
idx2 = 1092-1

m = -1.0 / (idx2-idx1)

q = -m * idx2

ii = np.arange(0,nx,1)

f = m * ii + q 
f = np.clip(f,a_min=0,a_max=1)

plt.figure()
plt.plot(ii,f,linestyle="solid")
plt.scatter(ii,f,s=1.5,zorder=0)
plt.show()




