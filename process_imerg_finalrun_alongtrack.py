# Process IMERG data
# 
# 30-minute data, 0.1 degree resolution
# combine IMERG files into a single dataset per year swapping dimensions so it is (time,lat,lon)
# take climatological mean over campaign time period in all years, write out
# save zonal mean over campaign time period in each year, write out

import xarray as xr
import numpy as np
import pandas as pd
import glob
import os
import time

start = time.time()

# time period of campaign
start_time = np.datetime64('2024-08-16T08:00:00')
end_time = np.datetime64('2024-09-23T22:59:00')

#Open DSHIP data
DSHIP = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/M203/Dship_data/data/meteor_meteo_dship_20240923.nc')

#Interpolate DSHIP data to hourly resolution
ship_time_interp = pd.date_range(start_time, end_time, freq='h')
ship_lat_interp = np.interp(ship_time_interp,DSHIP.time,DSHIP.lat)
ship_lon_interp = np.interp(ship_time_interp,DSHIP.time,DSHIP.lon)

# Filename base for IMERG data
filebase = "/huracan/tank4/cornell/ORCESTRA/imerg/final_run_V07_TropAtl_AugSep/3B-HHR.MS.MRG.3IMERG."

# For each year, extract data along Meteor's track
years = np.arange(2004,2017)  # 1998-2024
prec_alongtrack = np.full((len(years),len(ship_lat_interp)),np.nan)
iyear = 0
for yy in years:
    # print year
    print("Processing year:", yy)
    
    prec_alongtrack1 = np.full(len(ship_lat_interp), np.nan)  # Initialize array for one year
    
    # Load IMERG data
    file_paths = glob.glob(filebase + str(yy) + "*")
    
    #Filter out zero-size files
    file_paths = [f for f in file_paths if os.path.getsize(f) > 0]

    #Open the files and combine them into a single dataset
    IMERG_ds = xr.open_mfdataset(file_paths,combine='by_coords')
    
    end = time.time()
    print("Elapsed time for loading IMERG data:", end - start, "seconds")
    
    start = time.time()

    # Swap lat & lon dimensions, so it is (time,lat,lon)
    IMERG = IMERG_ds.transpose('time','lat','lon',...)
    
    # Convert cftime.DatetimeJulian to pandas.DatetimeIndex
    imergtimes = IMERG['time'].values
    imergtimes_converted = pd.to_datetime([t.strftime('%Y-%m-%d %H:%M:%S') for t in imergtimes])
    IMERG['time'] = imergtimes_converted
    end = time.time()
    print("Elapsed time for transpose and converted time:", end - start, "seconds")
    
    start = time.time()

    # Make lat array from -20 to 40 with 0.25 degree spacing
    lat_target = np.arange(-20, 40.25, 0.25)

    # Make lon array from -80 to 20 with 0.25 degree spacing
    lon_target = np.arange(-80, 20.25, 0.25)
    
    # Interpolate IMERG from 0.10 degree to the 0.25 degree target arrays and coarsen to hourly on the hour
    imerg_coarse = IMERG.interp(lat=lat_target,lon=lon_target,method='linear')
    end = time.time()
    print("Elapsed time for interpolation to 0.25 degree grid:", end - start, "seconds")
    
    start = time.time()
    imerg_coarse = imerg_coarse.coarsen(time=2,boundary='pad').mean()
    end = time.time()
    print("Elapsed time for coarsening to hourly:", end - start, "seconds")
    
    start = time.time()
    #Change year in ship track to this year
    ship_time_interp_newyear = ship_time_interp.map(lambda t: t.replace(year=yy))
    imerg_coarse = imerg_coarse.interp(lat=imerg_coarse.lat,lon=imerg_coarse.lon,time=ship_time_interp_newyear,method='linear')
    end = time.time()
    print("Elapsed time for changing year in ship track and interpolating to right times:", end - start, "seconds")
    
    start = time.time()
    imerg_coarse.load()  # Load the data into memory
    end = time.time()
    print("Elapsed time for loading IMERG data into memory:", end - start, "seconds")
    
    start = time.time()
    for itime in range(0,len(ship_lat_interp)):
        prec_alongtrack1[itime] = imerg_coarse.precipitation.sel(lat=ship_lat_interp[itime],lon=ship_lon_interp[itime], time=ship_time_interp_newyear[itime], method='nearest')
    end = time.time()
    print("Elapsed time for extracting precipitation along ship track:", end - start, "seconds")
    
    start = time.time()
    # Store the precipitation data for this year
    prec_alongtrack[iyear, :] = prec_alongtrack1
    end = time.time()
    print("Elapsed time for storing precipitation data for year:", end - start, "seconds")
    
    start = time.time()
    # Close the dataset to free up memory
    IMERG_ds.close()
    imerg_coarse.close()
    end = time.time()
    print("Elapsed time for closing datasets:", end - start, "seconds")
    
    iyear += 1 # advance year counter

start = time.time()    
#remove year from ship_time_interp
ship_time = ship_time_interp.strftime('%m-%d %H:%M:%S')

# Create xarray DataArray for precipitation along track
prec_alongtrack = xr.DataArray(prec_alongtrack, coords=[years, ship_time], dims=['year', 'time'])
prec_alongtrack.name = 'precipitation'
prec_alongtrack.attrs['units'] = 'mm/hr'
prec_alongtrack.attrs['description'] = 'hourly IMERG precipitation (coarsened to 0.25 deg) along ship track'

# Make xarray dataset
prec_alongtrack_ds = prec_alongtrack.to_dataset()
prec_alongtrack_ds.attrs['description'] = 'IMERG precipitation along ship track'
prec_alongtrack_ds.attrs['source'] = 'IMERG v07 Final Run'
prec_alongtrack_ds.attrs['history'] = 'Created 2025-06-23 by Allison Wing'

# Add other variables
prec_alongtrack_ds['ship_lat'] = ('time', ship_lat_interp)
prec_alongtrack_ds['ship_lon'] = ('time', ship_lon_interp)

# Save to netcdf file
prec_alongtrack_ds.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/prec_alongtrack.nc', mode='w')


end = time.time()
print("Elapsed time for writing out data:", end - start, "seconds")