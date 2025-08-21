# Process ERA-5 data
# 
# hourly data, 0.25 degree resolution
# extract along ship track for all years

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

# Filename base for ERA-5 CWV data
filebase_CWV = "/huracan/tank4/cornell/ORCESTRA/era5/total_column_water_vapour/"

# Filename base for ERA-5 CAPE data
filebase_CAPE = "/huracan/tank4/cornell/ORCESTRA/era5/convective_available_potential_energy/"

# For each year, extract data along Meteor's track
years = np.arange(1996,2025)  # 1998-2024
cwv_alongtrack = np.full((len(years),len(ship_lat_interp)),np.nan)
cape_alongtrack = np.full((len(years),len(ship_lat_interp)),np.nan)
iyear = 0
for yy in years:
    # print year
    print("Processing year:", yy)
    
    cwv_alongtrack1 = np.full(len(ship_lat_interp), np.nan)  # Initialize array for one year
    cape_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    
    # Load ERA-5 CWV data
    file_paths_CWV_08 = glob.glob(filebase_CWV + str(yy) + '08/' + "*.nc")
    file_paths_CWV_09 = glob.glob(filebase_CWV + str(yy) + '09/' + "*.nc")
    file_paths_CWV = file_paths_CWV_08 + file_paths_CWV_09
    era5_cwv = xr.open_mfdataset(file_paths_CWV,combine='by_coords')

    # Load ERA-5 CAPE data
    file_paths_CAPE_08 = glob.glob(filebase_CAPE + str(yy) + '08/' + "*.nc")
    file_paths_CAPE_09 = glob.glob(filebase_CAPE + str(yy) + '09/' + "*.nc")
    file_paths_CAPE = file_paths_CAPE_08 + file_paths_CAPE_09
    era5_cape = xr.open_mfdataset(file_paths_CAPE,combine='by_coords')

    end = time.time()
    print("Elapsed time for loading ERA-5 data:", end - start, "seconds")
    
    start = time.time()
    #Change year in ship track to this year
    ship_time_interp_newyear = ship_time_interp.map(lambda t: t.replace(year=yy))
    end = time.time()
    print("Elapsed time for changing year in ship track and interpolating to right times:", end - start, "seconds")
    
    start = time.time()
    era5_cwv.load()  # Load the data into memory
    era5_cape.load()
    end = time.time()
    print("Elapsed time for loading ERA-5 data into memory:", end - start, "seconds")

    start = time.time()
    for itime in range(0,len(ship_lat_interp)):
        cape_alongtrack1[itime] = era5_cape.cape.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime], valid_time=ship_time_interp_newyear[itime],method='nearest')
        cwv_alongtrack1[itime] = era5_cwv.tcwv.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime], valid_time=ship_time_interp_newyear[itime],method='nearest')
    end = time.time()
    print("Elapsed time for extracting ERA-5 along ship track:", end - start, "seconds")

    start = time.time()
    # Store the ERA-5 data for this year
    cape_alongtrack[iyear, :] = cape_alongtrack1
    cwv_alongtrack[iyear, :] = cwv_alongtrack1
    end = time.time()
    print("Elapsed time for storing ERA-5 data for year:", end - start, "seconds")

    start = time.time()
    # Close the dataset to free up memory
    era5_cwv.close()
    era5_cape.close()
    end = time.time()
    print("Elapsed time for closing datasets:", end - start, "seconds")
    
    iyear += 1 # advance year counter

start = time.time()    

# Create xarray DataArray for along track variables
cape_alongtrack = xr.DataArray(cape_alongtrack, coords=[years, ship_time_interp], dims=['year', 'time'])
cape_alongtrack.name = 'cape'
cape_alongtrack.attrs['units'] = 'J/kg'
cape_alongtrack.attrs['description'] = 'hourly ERA-5 CAPE along ship track'

cwv_alongtrack = xr.DataArray(cwv_alongtrack, coords=[years, ship_time_interp], dims=['year', 'time'])
cwv_alongtrack.name = 'cwv'
cwv_alongtrack.attrs['units'] = 'kg/m^2'
cwv_alongtrack.attrs['description'] = 'hourly ERA-5 CWV along ship track'

# Make xarray dataset
era5_alongtrack_ds = xr.Dataset()
era5_alongtrack_ds['cape'] = cape_alongtrack
era5_alongtrack_ds['cwv'] = cwv_alongtrack
era5_alongtrack_ds.attrs['description'] = 'ERA-5 data extracted along ship track'
era5_alongtrack_ds.attrs['source'] = 'ERA-5'
era5_alongtrack_ds.attrs['history'] = 'Created 2025-08-21 by Allison Wing'

# Add other variables
era5_alongtrack_ds['ship_lat'] = ('time', ship_lat_interp)
era5_alongtrack_ds['ship_lon'] = ('time', ship_lon_interp)

# Save to netcdf file
era5_alongtrack_ds.to_netcdf('/home/awing/orcestra/data/ERA5_alongtrack_allyears.nc', mode='w')


end = time.time()
print("Elapsed time for writing out data:", end - start, "seconds")