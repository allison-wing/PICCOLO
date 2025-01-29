# Process IMERG data
# 
# combine IMERG files into a single dataset swapping dimensions so it is (time,lat,lon), 
# restrict region, remove extra variables, write out reduced file

import xarray as xr
import glob
import os

# Load IMERG data
file_paths = glob.glob("/huracan/tank1/work/ORCESTRA/imerg/late_run_V07/3B-HHR-L.MS.MRG.3IMERG*")

#Open the files and combine them into a single dataset
IMERG_ds = xr.open_mfdataset(file_paths,engine='h5netcdf',group='Grid',combine='by_coords')

# Swap lat & lon dimensions, so it is (time,lat,lon)
IMERG_ds = IMERG_ds.transpose('time','lat','lon',...)

# Restrict the region
IMERG_ds = IMERG_ds.sel(lat=slice(-10,50),lon=slice(-90,60))

# Remove extra variables
IMERG_reduced = IMERG_ds.drop_vars(['probabilityLiquidPrecipitation','randomError','precipitationQualityIndex'])

#Write out to new netcdf file
IMERG_reduced.to_netcdf('/huracan/tank1/work/ORCESTRA/imerg/imerg_20240809.nc')
