import xarray as xr
import numpy as np
import pandas as pd
import glob
import os
import time

start = time.time()

yy = 2024  # year for which to process data

# time period of campaign
start_time = np.datetime64('2024-08-16T08:00:00')
end_time = np.datetime64('2024-09-23T22:59:00')

#Open DSHIP data
DSHIP = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/M203/Dship_data/data/meteor_meteo_dship_20240923.nc')

#Interpolate DSHIP data to hourly resolution
ship_time_interp = pd.date_range(start_time, end_time, freq='h')
ship_lat_interp = np.interp(ship_time_interp,DSHIP.time,DSHIP.lat)
ship_lon_interp = np.interp(ship_time_interp,DSHIP.time,DSHIP.lon)

# Make lat array from -20 to 40 with 0.25 degree spacing
lat_target = np.arange(-20, 40.25, 0.25)

# Make lon array from -80 to 20 with 0.25 degree spacing
lon_target = np.arange(-80, 20.25, 0.25)

# load imerg data
imerg = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/imerg/imerg_finalrun_' + str(yy) + '0809.nc')

# Convert cftime.DatetimeJulian to pandas.DatetimeIndex
imergtimes = imerg['time'].values
imergtimes_converted = pd.to_datetime([t.strftime('%Y-%m-%d %H:%M:%S') for t in imergtimes])
imerg['time'] = imergtimes_converted

# put imerg data on 0.25 degree target grid
imerg_grid = imerg.interp(lat=lat_target,lon=lon_target,method='linear') #interpolate to mimic grid
imerg_grid_coarse = imerg_grid.coarsen(time=2,boundary='pad').mean() #coarsen to hourly

# change time in ship track to this year
ship_time_interp_newyear = ship_time_interp.map(lambda t: t.replace(year=yy))
imerg_grid_coarse = imerg_grid_coarse.interp(lat=imerg_grid_coarse.lat,lon=imerg_grid_coarse.lon,time=ship_time_interp_newyear,method='linear') #interpolate to mimic times. For some reason I now (01/31/25) need to include lat and lon to keep them as indexes

# extract data along ship track
prec_alongtrack_yy = np.full(len(ship_lat_interp),np.nan)

for itime in range(0,len(ship_lat_interp)):
    prec_alongtrack_yy[itime] = imerg_grid_coarse.precipitation.sel(lat=ship_lat_interp[itime],lon=ship_lon_interp[itime], time=ship_time_interp_newyear[itime],method='nearest')

#remove year from ship_time_interp
ship_time = ship_time_interp.strftime('%m-%d %H:%M:%S')

# create xarray data array
prec_alongtrack_yy = xr.DataArray(prec_alongtrack_yy, coords=[ship_time_interp_newyear], dims=['time'])
prec_alongtrack_yy.name = 'precipitation'
prec_alongtrack_yy.attrs['units'] = 'mm/h'
prec_alongtrack_yy.attrs['description'] = 'IMERG precipitation along ship track for year ' + str(yy)

# make xarray dataset
prec_alongtrack_yy_ds = prec_alongtrack_yy.to_dataset()
prec_alongtrack_yy_ds.attrs['description'] = 'IMERG precipitation along ship track for year ' + str(yy)
prec_alongtrack_yy_ds.attrs['source'] = 'IMERG Final Run data'  

# add other variables
prec_alongtrack_yy_ds['ship_lat'] = ('time', ship_lat_interp)
prec_alongtrack_yy_ds['ship_lon'] = ('time', ship_lon_interp)

# save to netcdf file
prec_alongtrack_yy_ds.to_netcdf('/huracan/tank4/cornell/ORCESTRA/imerg/prec_alongtrack_' + str(yy) + '.nc')

end = time.time()
print("Elapsed time for processing IMERG data:", end - start, "seconds")