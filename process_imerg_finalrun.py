# Process IMERG data
# 
# combine IMERG files into a single dataset swapping dimensions so it is (time,lat,lon), 
# write out file for one year

import xarray as xr
import glob
import os
import time

yy = 2024

start = time.time()

# Load IMERG data (30-minute data, 0.1 degree resolution)
file_paths = glob.glob("/huracan/tank4/cornell/ORCESTRA/imerg/final_run_V07_TropAtl_AugSep/3B-HHR.MS.MRG.3IMERG." + str(yy) + "*")
#Filter out zero-size files
file_paths = [f for f in file_paths if os.path.getsize(f) > 0]

end = time.time()
print("Elapsed time for loading file paths:", end - start, "seconds")

start = time.time()
#Open the files and combine them into a single dataset
#IMERG_ds = xr.open_mfdataset(file_paths,engine='h5netcdf',group='Grid',combine='by_coords')
IMERG_ds = xr.open_mfdataset(file_paths,combine='by_coords')

end = time.time()
print("Elapsed time for opening dataset:", end - start, "seconds")

start = time.time()
# Swap lat & lon dimensions, so it is (time,lat,lon)
IMERG = IMERG_ds.transpose('time','lat','lon',...)

end = time.time()
print("Elapsed time for taking transpose:", end - start, "seconds")

IMERG_reduced = IMERG.drop_dims(['latv','lonv','nv'])

start = time.time()
#load data into memory before writing
IMERG_reduced.load()

end = time.time()
print("Elapsed time for loading data into memory:", end - start, "seconds")

start = time.time()
#Write out to new netcdf file
IMERG_reduced.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/imerg_finalrun_' + str(yy) + '0809.nc')

end = time.time()
print("Elapsed time for writing to file:", end - start, "seconds")
