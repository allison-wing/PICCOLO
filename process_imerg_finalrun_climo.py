# Process IMERG data
# 
# combine IMERG files into a single dataset swapping dimensions so it is (time,lat,lon), 
# write out file for one year

import xarray as xr
import numpy as np
import glob
import os

lonMin, lonMax = -62.1, -9.9
latMin, latMax = -2.1, 22.1

start_mon = 'Aug.'
start_day = '10'
start_time = '08-'+start_day+'T00:00:00'
end_mon = 'Sep.'
end_day = '30'
end_time = '09-'+end_day+'T00:00:00'

# Filename base
filebase = "/huracan/tank4/cornell/ORCESTRA/imerg/final_run_V07_TropAtl_AugSep/3B-HHR.MS.MRG.3IMERG."

# For each year, extract for time period of campaign and region
years = np.arange(1999,2022)  # 1999-2021
for yy in years:
    # Load IMERG data
    file_paths = glob.glob(filebase + str(yy) + "*")
    
    #Filter out zero-size files
    file_paths = [f for f in file_paths if os.path.getsize(f) > 0]

    #Open the files and combine them into a single dataset
    IMERG_ds = xr.open_mfdataset(file_paths,combine='by_coords')

    # Swap lat & lon dimensions, so it is (time,lat,lon)
    IMERG = IMERG_ds.transpose('time','lat','lon',...)
    
    # Select region and time period of interest
    IMERGyear = IMERG.sel(latitude=slice(latMin,latMax),longitude=slice(lonMin,lonMax),time=slice(str(yy)+'-'+start_time, str(yy)+'-'+end_time))
    if yy == 1999:
        IMERGCampaign = IMERGyear
    else:
        IMERGCampaign = xr.concat([IMERGCampaign,IMERGyear],dim='time')

#Take climatological mean over campaign time period over all years
IMERGClimoMean = IMERGCampaign.mean(dim='time')

#Save to netcdf file
IMERGClimoMean.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/imerg_campaign_climo.nc')