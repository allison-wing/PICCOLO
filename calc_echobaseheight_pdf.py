## Plot echo top height for each frame of SEA-POL volume data 
# Restrict to AP

import xarray as xr
import numpy as np
import pandas as pd

from scipy.interpolate import interp2d, RectBivariateSpline
from datetime import datetime, timedelta
import cftime
import json
import glob
import os

# read in data
seapol = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4/PICCOLO_level4_volume_3D.nc')

# set to nan outside of radius 120 km to only include data with the 3D volume
radius_outer = 120  # km
radius_inner = 50

#find the min height where the reflectivity is above a threshold
threshold = -10

#Define time period for spatial map
#APtime = np.datetime64('2024-08-28T20:20:00')
APtime = np.datetime64('2024-08-17T00:00:00')
indexAP = np.abs(pd.to_datetime(seapol.time) - APtime).argmin()

# Make regular 10-minute time series
end_time = np.datetime64('2024-09-23T16:50:00')
time10m = pd.date_range(APtime, end_time, freq='10 min')
time10m = pd.to_datetime(time10m)

# Loop over times since AP
#for i in range(0, np.size(seapol.time[indexAP:-1])): #loop over all times since AP
for i in range(0,len(time10m)):  # loop over all times in the 10-minute series since AP (VOL1 only)

    #check if there is data for this time
    #if seapol.time[i].values not in seapol.time.values:
    if time10m[i] not in seapol.time.values: #VOL1 only
        #print(f"No data for time {seapol.time[i].values}")
        print(f"No data for time {time10m[i]}") #VOL1 only
        continue
    else:
        #print('Time:', seapol.time[indexAP + i].values)
        #map_dbz = seapol.DBZ[indexAP+i,:,:,:]
         
        print('Time:', seapol.time.sel(time=time10m[i]).values) #VOL1 only
        map_dbz = seapol.DBZ.sel(time=time10m[i]) #VOL1 only

        # set to nan outside of radius_outer to only include data with the 3D volume
        distances = np.sqrt((map_dbz.latitude - map_dbz.latitude[120, 120])**2 + (map_dbz.longitude - map_dbz.longitude[120, 120])**2) * 111.32  # Approximate conversion from degrees to km
        map_dbz = map_dbz.where(distances<=radius_outer,np.nan)  # Set values outside the radius to NaN

        # also mask within radius_inner, to limit to radial range where we have "full" 3D coverage
        map_dbz = map_dbz.where(distances >= radius_inner, np.nan)

        # mask for valid (non-NaN) data
        valid_data = ~np.isnan(map_dbz.values)

        # mask for reflectivity above threshold
        above_thresh = valid_data & (map_dbz.values >= threshold)

        # Find the lowest index (height) where above_thresh is True for each (y, x)
        min_indices = np.where(above_thresh.any(axis=0), np.argmax(above_thresh, axis=0), -1)

        # Initialize output: NaN where no valid data, -5 where threshold not met
        echo_base_height = np.full(map_dbz.shape[1:], np.nan)
        has_valid = valid_data.any(axis=0)
        echo_base_height[has_valid] = -5

        # Set echo base height where threshold is met
        valid = min_indices != -1
        echo_base_height[valid] = map_dbz.Z.values[min_indices[valid].astype(int)]

        # Reshape valid echo base height into a 1D array
        echo_base_height_flat = echo_base_height[valid].flatten()

        # Concatenate the echo base heights into a single array
        if i == 0:
            all_echo_base_heights = echo_base_height_flat
        else:
            all_echo_base_heights = np.concatenate((all_echo_base_heights, echo_base_height_flat))
    

# Save the echo base heights to netcdf
ds = xr.Dataset({'echo_base_height': (['data points'], all_echo_base_heights)},
                coords={'data points': np.arange(len(all_echo_base_heights))})
ds.to_netcdf('../../data/SEA-POLv1.0_echo_base_height_vol1_50_120_-10dbz.nc')