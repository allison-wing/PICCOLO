# %% [markdown]
# # SEA-POL low-level gridded rain rate

# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, colors, ticker
import matplotlib.dates as mdates
from scipy.interpolate import interp2d, RectBivariateSpline
from datetime import datetime, timedelta
import pandas as pd
import cftime
import seaborn as sns
import json
from matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap,Normalize
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cmweather
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# %%
seapol = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4b/PICCOLO_level4b_rainrate_2D.nc')


# %%
# Mask out missing data (-32769 = no data possible)
rainrate = seapol.RAINRATE.where(seapol.RAINRATE >=-30000, np.nan)
dbz = seapol.DBZ.where(seapol.DBZ >=-30000, np.nan)

# %%
# Change -9999 missing data to zeros (data possible but removed = "not raining" though technically could be below beam) 
rainrate2 = rainrate.where(rainrate != -9999., 0)
#dbz2 = dbz
dbz2 = dbz.where(rainrate != -9999., np.nan)

# %%
# Also mask out seemingly bad data (rainrate > 10000 mm/h) or > 1000 mm/h --> set to zero (like saying there is no echo)
rainrate3 = rainrate2.where((rainrate2 <= 10000) | (rainrate2.isnull()), 0) # keep nans
dbz3 = dbz2.where((rainrate2 <= 10000) | (dbz2.isnull()), np.nan)
rainrate4 = rainrate3.where((rainrate3 <= 1000) | (rainrate3.isnull()), 0)
dbz4 = dbz3.where((rainrate3 <= 1000) | (dbz3.isnull()), np.nan)

# or donn't
#rainrate4 = rainrate2
#dbz4 = dbz2

# %% [markdown]
# # Take spatial averages

# %%
#Spatial mean over 245 km x 245 km
rain245 = rainrate4.mean(dim=('X','Y'),skipna=True)

#Spatial mean over 120 km x 120 km
rain120 = rainrate4.sel(X=slice(-120000,120000), Y=slice(-120000,120000)).mean(dim=('X','Y'),skipna=True)

#Spatial mean over 60 km x 60 km
rain60 = rainrate4.sel(X=slice(-60000,60000), Y=slice(-60000,60000)).mean(dim=('X','Y'),skipna=True)

#Spatial mean over 12 km x 12 km
rain12 = rainrate4.sel(X=slice(-12000,12000), Y=slice(-12000,12000)).mean(dim=('X','Y'),skipna=True)

#Spatial mean within 1 km
rain1 = rainrate4.sel(X=slice(-1000,1000), Y=slice(-1000,1000)).mean(dim=('X','Y'),skipna=True)


# %%
# Find indices of known bad data times (to mask out)
# Search times
time1a= np.datetime64('2024-08-16T00:00:00') #CV islands
time1b = np.datetime64('2024-08-17T02:00:00') #CV islands
time2a = np.datetime64('2024-08-27T22:00:00') #Praia
time2b = np.datetime64('2024-08-28T07:00:00') #Praia
time3 = np.datetime64('2024-09-21T18:10:00') #zeros are NaNs instead

#Find indices for start and end times
index1a = np.abs(pd.to_datetime(rain245.time) - time1a).argmin()
index1b = np.abs(pd.to_datetime(rain245.time) - time1b).argmin()
index2a = np.abs(pd.to_datetime(rain245.time) - time2a).argmin()
index2b = np.abs(pd.to_datetime(rain245.time) - time2b).argmin()
index3 = np.abs(pd.to_datetime(rain245.time) - time3).argmin()
print(f"Index for time 1a: {index1a}")
print(f"Index for time 1b: {index1b}")
print(f"Index for time 2a: {index2a}")
print(f"Index for time 2b: {index2b}")
print(f"Index for time 3: {index3}")


#Set bad data to NaN
rain245[index1a:index1b+1] = np.nan
rain245[index2a:index2b+1] = np.nan
rain245[index3] = np.nan
rain120[index1a:index1b+1] = np.nan
rain120[index2a:index2b+1] = np.nan
rain120[index3] = np.nan
rain60[index1a:index1b+1] = np.nan
rain60[index2a:index2b+1] = np.nan
rain60[index3] = np.nan
rain12[index1a:index1b+1] = np.nan
rain12[index2a:index2b+1] = np.nan
rain12[index3] = np.nan
rain1[index1a:index1b+1] = np.nan
rain1[index2a:index2b+1] = np.nan
rain1[index3] = np.nan



# %%
# Conditional mean (only where rainrate > 0)
rain245cond = rainrate4.where(rainrate4>0).mean(dim=('X','Y'),skipna=True)
rain120cond = rainrate4.sel(X=slice(-120000,120000), Y=slice(-120000,120000)).where(rainrate4.sel(X=slice(-120000,120000), Y=slice(-120000,120000))>0).mean(dim=('X','Y'),skipna=True)
rain60cond = rainrate4.sel(X=slice(-60000,60000), Y=slice(-60000,60000)).where(rainrate4.sel(X=slice(-60000,60000), Y=slice(-60000,60000))>0).mean(dim=('X','Y'),skipna=True)
rain12cond = rainrate4.sel(X=slice(-12000,12000), Y=slice(-12000,12000)).where(rainrate4.sel(X=slice(-12000,12000), Y=slice(-12000,12000))>0).mean(dim=('X','Y'),skipna=True)
rain1cond = rainrate4.sel(X=slice(-1000,1000), Y=slice(-1000,1000)).where(rainrate4.sel(X=slice(-1000,1000), Y=slice(-1000,1000))>0).mean(dim=('X','Y'),skipna=True)

# %%
#set bad data to NaN
rain245cond[index1a:index1b+1] = np.nan
rain245cond[index2a:index2b+1] = np.nan

rain120cond[index1a:index1b+1] = np.nan
rain120cond[index2a:index2b+1] = np.nan

rain60cond[index1a:index1b+1] = np.nan
rain60cond[index2a:index2b+1] = np.nan

rain12cond[index1a:index1b+1] = np.nan
rain12cond[index2a:index2b+1] = np.nan

rain1cond[index1a:index1b+1] = np.nan
rain1cond[index2a:index2b+1] = np.nan



# %% [markdown]
# # Calculate fractional area coverage of precip (>0)

# %%
# Find number of points with rainrate > 0
rain245count = rainrate4.where(rainrate4>0).count(dim=('X','Y'))
rain120count = rainrate4.sel(X=slice(-120000,120000), Y=slice(-120000,120000)).where(rainrate4.sel(X=slice(-120000,120000), Y=slice(-120000,120000))>0).count(dim=('X','Y'))
rain60count = rainrate4.sel(X=slice(-60000,60000), Y=slice(-60000,60000)).where(rainrate4.sel(X=slice(-60000,60000), Y=slice(-60000,60000))>0).count(dim=('X','Y'))
rain12count = rainrate4.sel(X=slice(-12000,12000), Y=slice(-12000,12000)).where(rainrate4.sel(X=slice(-12000,12000), Y=slice(-12000,12000))>0).count(dim=('X','Y'))
rain1count = rainrate4.sel(X=slice(-1000,1000), Y=slice(-1000,1000)).where(rainrate4.sel(X=slice(-1000,1000), Y=slice(-1000,1000))>0).count(dim=('X','Y'))  

# Calculate fractional area coverage of precip (>0)
rain245frac = rain245count / (rainrate4.X.size * rainrate4.Y.size)
rain120frac = rain120count / (rainrate4.sel(X=slice(-120000,120000), Y=slice(-120000,120000)).X.size * rainrate4.sel(X=slice(-120000,120000), Y=slice(-120000,120000)).Y.size)
rain60frac = rain60count / (rainrate4.sel(X=slice(-60000,60000), Y=slice(-60000,60000)).X.size * rainrate4.sel(X=slice(-60000,60000), Y=slice(-60000,60000)).Y.size)
rain12frac = rain12count / (rainrate4.sel(X=slice(-12000,12000), Y=slice(-12000,12000)).X.size * rainrate4.sel(X=slice(-12000,12000), Y=slice(-12000,12000)).Y.size)
rain1frac = rain1count / (rainrate4.sel(X=slice(-1000,1000), Y=slice(-1000,1000)).X.size * rainrate4.sel(X=slice(-1000,1000), Y=slice(-1000,1000)).Y.size)  

# Set bad data to NaN
rain245frac[index1a:index1b+1] = np.nan
rain120frac[index1a:index1b+1] = np.nan
rain60frac[index1a:index1b+1] = np.nan
rain12frac[index1a:index1b+1] = np.nan
rain1frac[index1a:index1b+1] = np.nan

# %% [markdown]
# ## Pad missing data with nans so it is on a regular 10-min time grid

# %%
# Make regular 10-minute time series
start_time = np.datetime64('2024-08-16T08:00:00')
end_time = np.datetime64('2024-09-23T16:50:00')
time10m = pd.date_range(start_time, end_time, freq='10 min')
time10m = pd.to_datetime(time10m)

# Make placeholder data array of nans of len(time10m)
rain245mean = np.full(len(time10m), np.nan)
rain120mean = np.full(len(time10m), np.nan)
rain60mean = np.full(len(time10m), np.nan)
rain12mean = np.full(len(time10m), np.nan)
rain1mean = np.full(len(time10m), np.nan)

rain245int = np.full(len(time10m), np.nan)
rain120int = np.full(len(time10m), np.nan)
rain60int = np.full(len(time10m), np.nan)
rain12int = np.full(len(time10m), np.nan)
rain1int = np.full(len(time10m), np.nan)

rain245fa = np.full(len(time10m), np.nan)
rain120fa = np.full(len(time10m), np.nan)
rain60fa = np.full(len(time10m), np.nan)
rain12fa = np.full(len(time10m), np.nan)
rain1fa = np.full(len(time10m), np.nan)

# Convert both time arrays to pandas.DatetimeIndex for reliable comparison
seapol_times = pd.to_datetime(seapol.time.values)

# Find indices in time10m that are exactly in seapol.time
mask = np.isin(time10m, seapol_times)
matched_idx = np.where(mask)[0]

# Find corresponding indices in seapol.time for the matched times
seapol_idx = {t: i for i, t in enumerate(seapol_times)}
matched_seapol_idx = [seapol_idx[t] for t in time10m[mask]]

# Fill arrays only at matched indices
rain245mean[matched_idx]    = rain245.values[matched_seapol_idx]
rain120mean[matched_idx]    = rain120.values[matched_seapol_idx]
rain60mean[matched_idx]     = rain60.values[matched_seapol_idx]
rain12mean[matched_idx]     = rain12.values[matched_seapol_idx]
rain1mean[matched_idx]      = rain1.values[matched_seapol_idx]

rain245int[matched_idx] = rain245cond.values[matched_seapol_idx]
rain120int[matched_idx] = rain120cond.values[matched_seapol_idx]
rain60int[matched_idx]  = rain60cond.values[matched_seapol_idx]
rain12int[matched_idx]  = rain12cond.values[matched_seapol_idx]
rain1int[matched_idx]   = rain1cond.values[matched_seapol_idx]

rain245fa[matched_idx] = rain245frac.values[matched_seapol_idx]
rain120fa[matched_idx] = rain120frac.values[matched_seapol_idx]
rain60fa[matched_idx]  = rain60frac.values[matched_seapol_idx]
rain12fa[matched_idx]  = rain12frac.values[matched_seapol_idx]
rain1fa[matched_idx]   = rain1frac.values[matched_seapol_idx]

# Convert to xarray data arrays
rain245mean = xr.DataArray(rain245mean, coords=[time10m], dims=['time'], name='rainrate_245km_mean')
rain120mean = xr.DataArray(rain120mean, coords=[time10m], dims=['time'], name='rainrate_120km_mean')
rain60mean = xr.DataArray(rain60mean, coords=[time10m], dims=['time'], name='rainrate_60km_mean')
rain12mean = xr.DataArray(rain12mean, coords=[time10m], dims=['time'], name='rainrate_12km_mean')
rain1mean = xr.DataArray(rain1mean, coords=[time10m], dims=['time'], name='rainrate_1km_mean')

rain245int = xr.DataArray(rain245int, coords=[time10m], dims=['time'], name='rainrate_245km_int')
rain120int = xr.DataArray(rain120int, coords=[time10m], dims=['time'], name='rainrate_120km_int')
rain60int = xr.DataArray(rain60int, coords=[time10m], dims=['time'], name='rainrate_60km_int')
rain12int = xr.DataArray(rain12int, coords=[time10m], dims=['time'], name='rainrate_12km_int')
rain1int = xr.DataArray(rain1int, coords=[time10m], dims=['time'], name='rainrate_1km_int')

rain245fa = xr.DataArray(rain245fa, coords=[time10m], dims=['time'], name='rainrate_245km_frac')
rain120fa = xr.DataArray(rain120fa, coords=[time10m], dims=['time'], name='rainrate_120km_frac')
rain60fa = xr.DataArray(rain60fa, coords=[time10m], dims=['time'], name='rainrate_60km_frac')
rain12fa = xr.DataArray(rain12fa, coords=[time10m], dims=['time'], name='rainrate_12km_frac')
rain1fa = xr.DataArray(rain1fa, coords=[time10m], dims=['time'], name='rainrate_1km_frac')


# %% [markdown]
# # Write out spatial averages to file

# %%
#combine into one dataset
#regular 10-min time array with nans where missing data
rainrate = xr.Dataset({'rain245_mean': rain245mean, 'rain120_mean': rain120mean, 'rain60_mean': rain60mean, 'rain12_mean': rain12mean, 'rain1_mean': rain1mean, 
                       'rain245_int': rain245int, 'rain120_int': rain120int, 'rain60_int': rain60int, 'rain12_int': rain12int, 'rain1_int': rain1int,
                       'rain245_frac': rain245fa, 'rain120_frac': rain120fa, 'rain60_frac': rain60fa, 'rain12_frac': rain12fa, 'rain1_frac': rain1fa})

#native time array that is blnk where missing data
#rainrate = xr.Dataset({'rain245_mean': rain245, 'rain120_mean': rain120, 'rain60_mean': rain60, 'rain12_mean': rain12, 'rain1_mean': rain1, 'rain245_int': rain245cond, 'rain120_int': rain120cond, 'rain60_int': rain60cond, 'rain12_int': rain12cond, 'rain1_int': rain1cond})

#add attributes
rainrate.attrs['title'] = 'Spatial mean rainrate from SEA-POL long-range, low-elevation scans'
rainrate.attrs['description'] = ('Spatial means of rainrate from SEA-POL long-range, low-elevation scans at different spatial scales. '
                                 'Mean is the average over all valid data points (including zeros). '
                                 'Int (intensity) is the average over all data points where rainrate > 0.'
                                 'fa is the fractional area coverage of rainrate > 0.')
rainrate.attrs['source'] = 'SEA-POL Level4b Gridded 2D Rain Rate'
rainrate.attrs['units'] = 'mm/h'

#save to netcdf
rainrate.to_netcdf('/home/awing/orcestra/data/SEA-POL_4b_rainrate_2D_spatial_means_masked_reg10.nc')