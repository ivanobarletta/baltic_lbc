import pandas as pd
from glob import glob 

# this script concatenates the NEMO volume_transport files
# that are present in the subfolders of rootDir
#
# rootDir/
#   R20211229/volume_transport
#   R20220105/volume_transport
#   R20220112/volume_transport
#   R20220119/volume_transport
#   R20220126/volume_transport

Sv  = 1e6   # m3/s conversion factor from Sverdrup (as in NEMO volume_transport)

# root directory of simulation
rootDir     = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/EXP_BT0/"

# initial day
date0       = "2021-12-29"
# number of weeks
NPeriods    = 3

sectionName = "KATTEGAT"

sectionName = "EAST2"
sectionName = "SKAGERRAK"
sectionName = "GIBRALTAR"

outFile     = "volume_transport_%s_BT0_diadct" % sectionName

# create a date range of weeks
dates   = pd.date_range(start=date0,periods=NPeriods,freq="7d")

# build list of folders
folders = [ rootDir + date.strftime("R%Y%m%d") for date in dates ] 

files   = [ folder+"/volume_transport" for folder in folders]

files   = glob("%s/*/volume_transport" % rootDir)

def read_transport_df(path="",sectionName="transectName"):

    # read dataframe     
    df = pd.read_csv(path,header=None,sep="\s+")

    print(df.info)
    print()

    # assign names to columns
    df.columns = ["date", "time_step", "section_number", "section_name", "section_slope_coefficient", "class_number", "class_name", "class_bound1" , 
                    "class_bound2", "transport_direction1" , "transport_direction2", "transport_total"]

    # filter by section_name and class_name (always total)
    df = df.loc[ (df["section_name"] == sectionName ) & ( df["class_name"] == "total" ) ]

    Nperiods = df.shape[0]

    # create dates 
    new_dates = pd.date_range(start=pd.to_datetime(df["date"].iloc[0],format="%Y%m%d"), periods=Nperiods,freq="1h") + pd.Timedelta("1h")

    # format new_dates
    new_dates = new_dates.strftime("%Y%m%dT%H%M%S")

    # replace dates in dataframe with new_ones
    df.index = new_dates

    # select only date,section_name,transport_total
    df  = df[["section_name","transport_total"]]

    # convert transport to m3/s
    df["transport_total"] = df["transport_total"].apply(lambda x: x*Sv)


    return df

"""
with open(outFile, 'w') as outStream:
    for fname in files:
        with open(fname) as infile:
            for line in infile:
                outStream.write(line)
"""

# concatenate all dataframes
print("Processing Files:")
for file in files:
    print("%s" % file)

dfOut   = pd.concat([read_transport_df(path=file,sectionName=sectionName) for file in files ])

# save to file
dfOut.to_csv(outFile,index=True)

# tips for reading the output file
# Do:
#   df = pd.read_csv("file",sep=",",index_col=0) 
#
#   this will consider the 1st colums as index    
#
# Then convert to pandas Series:
#
#   time_series = pd.Series(df["transport_total"].values, index=pd.to_datetime(df["date"],format="%Y%m%dT%H%M%S"))
#
# with the timeseries you can do, for example, daily average:
#
#   time_series.resample(rule=pd.Timedelta("1d")).mean().shift(0.5,freq="1d").plot(color="0.6",linestyle="dashdot")
#
# the half-day shift is necessary because the resampling refers the output values
# to the left of the bin rather to the center



