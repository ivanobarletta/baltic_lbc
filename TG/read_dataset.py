import pandas as pd

path    = "in_situ/index_monthly.txt"

varName = "SLEV"

# lon_min,lon_max,lat_min,lat_max
box     = [9,14,53,60]

df      = pd.read_csv(path,header=5)

print("initial length: %d"  % len(df))

print(df.columns)

# convert dates from string to datetime
df["time_coverage_start"]   = pd.to_datetime(df["time_coverage_start"],format="%Y-%m-%dT%H:%M:%SZ")
df["time_coverage_end"]     = pd.to_datetime(df["time_coverage_end"],format="%Y-%m-%dT%H:%M:%SZ")

# filter by dates 
df  = df.loc[ (df["time_coverage_start"] > "2022-01-01" ) & (df["time_coverage_end"] < "2023-12-31")]

# filter by variable name
df  = df[df["parameters"].str.contains(varName)]

# filter by geographical area
df  = df[(df["geospatial_lon_min"] > box[0]) & (df["geospatial_lon_max"] < box[1]) & (df["geospatial_lat_min"] > box[2]) & (df["geospatial_lat_max"] < box[3])]

print("length of resulting df: %d" % len(df))
