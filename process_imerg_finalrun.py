# Process IMERG data
# 
# combine IMERG files into a single dataset swapping dimensions so it is (time,lat,lon), 
# write out file for one year

import xarray as xr
import glob
import os

# Load IMERG data
file_paths = glob.glob("/huracan/tank4/cornell/ORCESTRA/imerg/final_run_V07_TropAtl_AugSep/3B-HHR.MS.MRG.3IMERG.2024*")

#Open the files and combine them into a single dataset
#IMERG_ds = xr.open_mfdataset(file_paths,engine='h5netcdf',group='Grid',combine='by_coords')
IMERG_ds = xr.open_mfdataset(file_paths,combine='by_coords')

# Swap lat & lon dimensions, so it is (time,lat,lon)
IMERG = IMERG_ds.transpose('time','lat','lon',...)

#Write out to new netcdf file
IMERG.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/imerg_finalrun_20240809.nc')