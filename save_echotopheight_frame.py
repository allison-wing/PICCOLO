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

import matplotlib.pyplot as plt
from matplotlib import rc, colors, ticker
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import seaborn as sns

# Load colormap

# Read the colormap from a .txt file
def load_colormap_from_txt(file_path):
    # Load RGB values from the file
    rgb_values = np.loadtxt(file_path)
    return ListedColormap(rgb_values)

# Example usage
colormap_file = 'chase-spectral-rgb.txt'  # Replace with your .txt file path
radar_cmap = load_colormap_from_txt(colormap_file)

discrete_cmap = ListedColormap(radar_cmap(np.linspace(0, 1, 16)))
discrete_mask_cmap = discrete_cmap.copy()
discrete_mask_cmap.set_under(color='white') #set values below vmin to white
discrete_mask_cmap.set_bad(color='lightgray') #set missing (NaN) values to gray

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
    print('Processing frame:', i)
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
    
    # plot
    fig, axs = plt.subplots(figsize=(8, 8),layout="constrained")  
    cax = axs.pcolormesh(map_dbz.longitude, map_dbz.latitude, echo_top_height/1000, cmap=discrete_mask_cmap, vmin=0, vmax=15)
    axs.set_aspect('equal',adjustable='box')
    cbar = fig.colorbar(cax, ax=axs, orientation='vertical',pad=0.02, shrink=0.70)
    cbar.ax.tick_params(labelsize=16)
    #cbar.set_ticks(np.linspace(-10,60,8))
    cbar.ax.set_ylabel('Echo Top Height (km)', fontsize=16)
    axs.set_title('SEA-POL: ' + str(seapol.time[indexAP + i].values),fontsize=16)
    axs.tick_params(axis='both', labelsize=16)
    axs.set_xlabel('Longitude', fontsize=16)
    axs.set_ylabel('Latitude', fontsize=16)

    # Add grid lines 
    axs.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    # Add 120 km circle
    axs.add_patch(plt.Circle((map_dbz.longitude[120,120], map_dbz.latitude[120,120]), 120/111.32, color='gray', alpha = 0.7,fill=False, linestyle='--', linewidth=0.5))
    # Add 60 km circle
    axs.add_patch(plt.Circle((map_dbz.longitude[120,120], map_dbz.latitude[120,120]), 60/111.32, color='gray', alpha = 0.7,fill=False, linestyle='--', linewidth=0.5))
    # Add 18 km circle
    axs.add_patch(plt.Circle((map_dbz.longitude[120,120], map_dbz.latitude[120,120]), 18/111.32, color='gray', alpha = 0.7,fill=False, linestyle='--', linewidth=0.5))
    
    #save figure
    fig.savefig(f'/home/awing/orcestra/figures/echo_top_height_frames/echo_top_height_{i:04d}.png', dpi=300, bbox_inches='tight')

