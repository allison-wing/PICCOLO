# Process IMERG data
# 
# 30-minute data, 0.1 degree resolution
# combine IMERG files into a single dataset per year swapping dimensions so it is (time,lat,lon)
# take climatological mean over campaign time period in all years, write out
# save zonal mean over campaign time period in each year, write out

import xarray as xr
import numpy as np
import glob
import os
import time

start = time.time()

# region of interest
lonMin, lonMax = -62.1, -9.9
latMin, latMax = -2.1, 22.1

# for zonal mean
lon1 = -30
lon2 = -20

# time period of campaign
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
    IMERGyear = IMERG.sel(lat=slice(latMin,latMax),lon=slice(lonMin,lonMax),time=slice(str(yy)+'-'+start_time, str(yy)+'-'+end_time)).mean(dim='time')
    zonalmean_precip_year = IMERGyear.sel(lon=slice(lon1,lon2)).mean(dim='lon')
    if yy == 1999:
        IMERGCampaign = IMERGyear
        zonalmean_precip = zonalmean_precip_year
    else:
        IMERGCampaign = xr.concat([IMERGCampaign,IMERGyear],dim='time')
        zonalmean_precip = xr.concat([zonalmean_precip,zonalmean_precip_year],dim='time')

#Take climatological mean over campaign time period over all years
IMERGClimoMean = IMERGCampaign.mean(dim='time')

#Save to netcdf file (climatological mean at each grid point over campaign time period in all years)
IMERGClimoMean.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/imerg_campaign_climo.nc')

# Save zonal mean to netcdf file (zonal mean over campaign time period in each year)
zonalmean_precip = zonalmean_precip.assign_coords(time=years)
zonalmean_precip.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/imerg_zonalmean_precip.nc')

end = time.time()
print("Elapsed time for processing IMERG data:", end - start, "seconds")