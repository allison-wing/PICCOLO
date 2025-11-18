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

# time period of campaign [really 8 UTC on 08-16 and 22:59 on 09-23, but adjusted here for era5 6-hourly times]
start_time = np.datetime64('2024-08-16T06:00:00')
end_time = np.datetime64('2024-09-24T00:00:00')

#Open DSHIP data
DSHIP = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/M203/Dship_data/data/meteor_meteo_dship_20240923.nc')

#Interpolate DSHIP data to 6-hourly resolution
ship_time_interp = pd.date_range(start_time, end_time, freq='6h')
ship_lat_interp = np.interp(ship_time_interp,DSHIP.time,DSHIP.lat)
ship_lon_interp = np.interp(ship_time_interp,DSHIP.time,DSHIP.lon)

# Filename base for ERA-5 10m wind data
filebase_10m = "/mars/tank3/era5/uv_10m_wind/"

start = time.time()
# Load ERA-5 u data because all years are in one file
era5_u850 = xr.open_dataset('/mars/tank3/era5/pressure_variables/u_component_of_wind/u_850.nc')
era5_u600 = xr.open_dataset('/mars/tank3/era5/pressure_variables/u_component_of_wind/u_600.nc')

# Load ERA-5 v data because all years are in one file
era5_v850 = xr.open_dataset('/mars/tank3/era5/pressure_variables/v_component_of_wind/v_850.nc')
era5_v600 = xr.open_dataset('/mars/tank3/era5/pressure_variables/v_component_of_wind/v_600.nc')
end = time.time()
print("Elapsed time for loading ERA-5 data:", end - start, "seconds")

# For each year, extract data along Meteor's track
years = np.arange(1990,2019) 
shear06_alongtrack = np.full((len(years),len(ship_lat_interp)),np.nan)
shear02_alongtrack = np.full((len(years),len(ship_lat_interp)),np.nan)
iyear = 0
for yy in years:
    # print year
    print("Processing year:", yy)
    
    shear06_alongtrack1 = np.full(len(ship_lat_interp), np.nan)  # Initialize array for one year
    shear02_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    u600_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    v600_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    u850_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    v850_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    u10m_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    v10m_alongtrack1 = np.full(len(ship_lat_interp),np.nan)
    
    # Load ERA-5 10m wind data
    file_paths_10m_08 = glob.glob(filebase_10m + str(yy) + '-08.nc')
    file_paths_10m_09 = glob.glob(filebase_10m + str(yy) + '-09.nc')
    file_paths_10m = file_paths_10m_08 + file_paths_10m_09
    era5_10m = xr.open_mfdataset(file_paths_10m,combine='by_coords')

    end = time.time()
    print("Elapsed time for loading ERA-5 data:", end - start, "seconds")
    
    start = time.time()
    #Change year in ship track to this year
    ship_time_interp_newyear = ship_time_interp.map(lambda t: t.replace(year=yy))
    end = time.time()
    print("Elapsed time for changing year in ship track:", end - start, "seconds")
    
    start = time.time()
    era5_10m.load()  # Load the data into memory
    end = time.time()
    print("Elapsed time for loading ERA-5 data into memory:", end - start, "seconds")

    start = time.time()
    for itime in range(0,len(ship_lat_interp)):
        u600_alongtrack1[itime] = era5_u600.u.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime]+360, time=ship_time_interp_newyear[itime],method='nearest')
        v600_alongtrack1[itime] = era5_v600.v.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime]+360, time=ship_time_interp_newyear[itime],method='nearest')
        u850_alongtrack1[itime] = era5_u850.u.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime]+360, time=ship_time_interp_newyear[itime],method='nearest')
        v850_alongtrack1[itime] = era5_v850.v.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime]+360, time=ship_time_interp_newyear[itime],method='nearest')
        u10m_alongtrack1[itime] = era5_10m.u10.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime]+360, time=ship_time_interp_newyear[itime],method='nearest')
        v10m_alongtrack1[itime] = era5_10m.v10.sel(latitude=ship_lat_interp[itime],longitude=ship_lon_interp[itime]+360, time=ship_time_interp_newyear[itime],method='nearest')
    end = time.time()
    print("Elapsed time for extracting ERA-5 along ship track:", end - start, "seconds")
    
    start = time.time()
    # calculate wind shear
    # surface to 600 hPa [approx equivalent of 0-6 km]
    shear06_alongtrack1 = np.sqrt( (u600_alongtrack1 - u10m_alongtrack1)**2 + (v600_alongtrack1 - v10m_alongtrack1)**2 )

    # surface to 850 hPa [approx equivalent of 0-2 km]
    shear02_alongtrack1 = np.sqrt( (u850_alongtrack1 - u10m_alongtrack1)**2 + (v850_alongtrack1 - v10m_alongtrack1)**2 )
    end = time.time()
    print("Elapsed time for calculating wind shear:", end - start, "seconds")

    start = time.time()
    # Store the ERA-5 data for this year
    shear06_alongtrack[iyear, :] = shear06_alongtrack1
    shear02_alongtrack[iyear, :] = shear02_alongtrack1
    end = time.time()
    print("Elapsed time for storing ERA-5 data for year:", end - start, "seconds")

    start = time.time()
    # Close the dataset to free up memory
    era5_10m.close()
    end = time.time()
    print("Elapsed time for closing datasets:", end - start, "seconds")
    
    iyear += 1 # advance year counter

start = time.time()
# Close the dataset to free up memory
era5_u600.close()
era5_v600.close()
era5_u850.close()
era5_v850.close()
end = time.time()
print("Elapsed time for closing datasets:", end - start, "seconds")

start = time.time()    

# Create xarray DataArray for along track variables
shear06_alongtrack = xr.DataArray(shear06_alongtrack, coords=[years, ship_time_interp], dims=['year', 'time'])
shear06_alongtrack.name = 'Shear: Sfc to 600 hPa'
shear06_alongtrack.attrs['units'] = 'm/s'
shear06_alongtrack.attrs['description'] = 'Deep shear along ship track'

shear02_alongtrack = xr.DataArray(shear02_alongtrack, coords=[years, ship_time_interp], dims=['year', 'time'])
shear02_alongtrack.name = 'Shear: Sfc to 850 hPa'
shear02_alongtrack.attrs['units'] = 'm/s'
shear02_alongtrack.attrs['description'] = 'Low-level shear along ship track'

# Make xarray dataset
era5_alongtrack_ds = xr.Dataset()
era5_alongtrack_ds['shear06'] = shear06_alongtrack
era5_alongtrack_ds['shear02'] = shear02_alongtrack
era5_alongtrack_ds.attrs['description'] = 'ERA-5 data extracted along ship track'
era5_alongtrack_ds.attrs['source'] = 'ERA-5'
era5_alongtrack_ds.attrs['history'] = 'Created 2025-11-17 by Allison Wing'

# Add other variables
era5_alongtrack_ds['ship_lat'] = ('time', ship_lat_interp)
era5_alongtrack_ds['ship_lon'] = ('time', ship_lon_interp)

# Save to netcdf file
era5_alongtrack_ds.to_netcdf('/home/awing/orcestra/data/ERA5_shear_alongtrack_allyears.nc', mode='w')


end = time.time()
print("Elapsed time for writing out data:", end - start, "seconds")