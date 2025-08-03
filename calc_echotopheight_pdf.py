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
seapol = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4b/PICCOLO_level4b_volume_3D.nc')

# set to nan outside of radius 120 km to only include data with the 3D volume
radius = 120  # km

#find the maximum height where the reflectivity is above a threshold
threshold = 10

#Define time period for spatial map
APtime = np.datetime64('2024-08-28T20:00:00')
indexAP = np.abs(pd.to_datetime(seapol.time) - APtime).argmin()

# Loop over times since AP
for i in range(0, np.size(seapol.time[indexAP:-1]) + 1):
    print('Time:', seapol.time[indexAP + i].values)
    map_dbz = seapol.DBZ[indexAP+i,:,:,:]

    # set to nan outside of radius 120 km to only include data with the 3D volume
    distances = np.sqrt((map_dbz.latitude - map_dbz.latitude[120, 120])**2 + (map_dbz.longitude - map_dbz.longitude[120, 120])**2) * 111.32  # Approximate conversion from degrees to km
    map_dbz = map_dbz.where(distances<=radius,np.nan)  # Set values outside the radius to NaN

    # mask for valid (non-NaN) data
    valid_data = ~np.isnan(map_dbz.values)

    # mask for reflectivity above threshold
    above_thresh = valid_data & (map_dbz.values >= threshold)

    # Find the highest index (height) where above_thresh is True for each (y, x)
    max_indices = np.where(above_thresh.any(axis=0), above_thresh.shape[0] - 1 - np.argmax(above_thresh[::-1], axis=0), -1)

    # Initialize output: NaN where no valid data, -5 where threshold not met
    echo_top_height = np.full(map_dbz.shape[1:], np.nan)
    has_valid = valid_data.any(axis=0)
    echo_top_height[has_valid] = -5

    # Set echo top height where threshold is met
    valid = max_indices != -1
    echo_top_height[valid] = map_dbz.Z.values[max_indices[valid].astype(int)]
    
    # Reshape valid echo top height into a 1D array
    echo_top_height_flat = echo_top_height[valid].flatten()
    
    # Concatenate the echo top heights into a single array
    if i == 0:
        all_echo_top_heights = echo_top_height_flat
    else:
        all_echo_top_heights = np.concatenate((all_echo_top_heights, echo_top_height_flat))
    

# Save the echo top heights to netcdf
ds = xr.Dataset({'echo_top_height': (['data points'], all_echo_top_heights)},
                coords={'data points': np.arange(len(all_echo_top_heights))})
ds.to_netcdf('../../data/SEA-POL_echo_top_height.nc')