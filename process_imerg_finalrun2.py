# Process IMERG data
# 
# combine IMERG files into a single dataset swapping dimensions so it is (time,lat,lon), 
# write out file for one year

import xarray as xr
import glob
import os
import time

yy = 2003

start = time.time()

# Load IMERG data (30-minute data, 0.1 degree resolution)
file_paths = glob.glob("/huracan/tank4/cornell/ORCESTRA/imerg/final_run_V07_TropAtl_AugSep/3B-HHR.MS.MRG.3IMERG." + str(yy) + "*")

#Open the files and combine them into a single dataset
#IMERG_ds = xr.open_mfdataset(file_paths,engine='h5netcdf',group='Grid',combine='by_coords')
IMERG_ds = xr.open_mfdataset(file_paths,combine='by_coords')

# Swap lat & lon dimensions, so it is (time,lat,lon)
IMERG = IMERG_ds.transpose('time','lat','lon',...)

#Write out to new netcdf file
IMERG.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/imerg_finalrun_' + str(yy) + '.nc')

end = time.time()
print("Elapsed time for processing IMERG data:", end - start, "seconds")
