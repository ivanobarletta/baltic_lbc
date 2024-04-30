import pandas as pd
import subprocess
import sys

#Purpose of the script:
#    1) read a .csv inventory that contains a list of locations
#    2) run the script (defined in sbatchFile ) to estract a timeseries for each one of
#        the locations


if len(sys.argv) < 2:
    raise Exception("Error: Provide path of inventory (csv")
    sys.exit(1)

pathInventory = sys.argv[1]

sbatchFile = "extract_ssh_ncks2.sh"

df = pd.read_csv(pathInventory,header=0,sep=",")

for fileName,lon,lat in zip(df["filename"],df["lon"],df["lat"]):
    fileName2 = fileName.replace(".nc","")
    print (lon,lat,fileName2)
    out_subset = subprocess.run(["sbatch",sbatchFile,"%s" % lat,"%s"%lon,"%s" % fileName2])
