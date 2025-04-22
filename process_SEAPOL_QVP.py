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
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cmweather
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import glob
import os

################ Before Praia ################
# Load data
file_dirs = glob.glob("/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4a/qvp_volume/202408*")

# loop through directories
for file_dir in file_dirs:
    # get the date from the directory name
    date = os.path.basename(file_dir)
    print(f"Processing data for {date}")

    file_paths = glob.glob(file_dir+"/*.nc")
    qvp = xr.open_mfdataset(file_paths, combine='by_coords')

    # Mask the data
    # missing data (blocking sectors)
    #rainrate = qvp.RATE_CSU_BLENDED.where(qvp.RATE_CSU_BLENDED >=-30000, np.nan)
    dbz = qvp.DBZ.where(qvp.DBZ >=-30000, np.nan)

    #change -9999 missing data to zeros for rain, NaN for DBZ
    #rainrate2 = rainrate.where(rainrate != -9999, 0)
    dbz2 = dbz.where(dbz != -9999, np.nan)

    # Average spatially at each level (results in time x height)
    #rainrate_avg = rainrate2.mean(dim=['X', 'Y'],skipna=True)
    dbz_avg = dbz2.mean(dim=['X', 'Y'],skipna=True) 
    #rainrate_cond = rainrate2.where(rainrate2>0).mean(dim=('X','Y'),skipna=True)

    # Save to netcdf file
    #rainrate_avg.attrs['units'] = 'mm/h'
    #rainrate_cond.attrs['units'] = 'mm/h'
    dbz_avg.attrs['units'] = 'dBZ'

    #combine into one dataset
    #qvp1D= xr.Dataset({'rainrate_avg': rainrate_avg, 'dbz_avg': dbz_avg, 'rainrate_int': rainrate_cond})
    qvp1D= xr.Dataset({'dbz_avg': dbz_avg})

    #add attributes
    qvp1D.attrs['title'] = 'Average rain rate and reflectivity at each altitude from SEA-POL QVP'
    qvp1D.attrs['description'] = 'Spatial means of rainrate from SEA-POL 45 degree elevation scans at different altitudes. Mean is the average over all valid data points (including zeros). Int (intensity) is the average over all data points where rainrate > 0.'
    qvp1D.attrs['source'] = 'SEA-POL Level4a QVP Volume'

    #save to netcdf
    qvp1D.to_netcdf('/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4a/qvp_1D/qvp_1D_'+date+'.nc', mode='w', format='NETCDF4')

############## After Praia - September #################
# Load data
file_dirs = glob.glob("/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4a/qvp_volume/202409*")

# loop through directories
for file_dir in file_dirs:
    # get the date from the directory name
    date = os.path.basename(file_dir)
    print(f"Processing data for {date}")

    file_paths = glob.glob(file_dir+"/*.nc")
    qvp = xr.open_mfdataset(file_paths, combine='by_coords')

    # Mask the data
    # missing data (blocking sectors)
    rainrate = qvp.RATE_CSU_BLENDED.where(qvp.RATE_CSU_BLENDED >=-30000, np.nan)
    dbz = qvp.DBZ.where(qvp.DBZ >=-30000, np.nan)

    #change -9999 missing data to zeros for rain, NaN for DBZ
    rainrate2 = rainrate.where(rainrate != -9999, 0)
    dbz2 = dbz.where(dbz != -9999, np.nan)

    # Average spatially at each level (results in time x height)
    rainrate_avg = rainrate2.mean(dim=['X', 'Y'],skipna=True)
    dbz_avg = dbz2.mean(dim=['X', 'Y'],skipna=True) 
    rainrate_cond = rainrate2.where(rainrate2>0).mean(dim=('X','Y'),skipna=True)

    # Save to netcdf file
    rainrate_avg.attrs['units'] = 'mm/h'
    rainrate_cond.attrs['units'] = 'mm/h'
    dbz_avg.attrs['units'] = 'dBZ'

    #combine into one dataset
    qvp1D= xr.Dataset({'rainrate_avg': rainrate_avg, 'dbz_avg': dbz_avg, 'rainrate_int': rainrate_cond})
    #qvp1D= xr.Dataset({'dbz_avg': dbz_avg})

    #add attributes
    qvp1D.attrs['title'] = 'Average rain rate and reflectivity at each altitude from SEA-POL QVP'
    qvp1D.attrs['description'] = 'Spatial means of rainrate from SEA-POL 45 degree elevation scans at different altitudes. Mean is the average over all valid data points (including zeros). Int (intensity) is the average over all data points where rainrate > 0.'
    qvp1D.attrs['source'] = 'SEA-POL Level4a QVP Volume'

    #save to netcdf
    qvp1D.to_netcdf('/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4a/qvp_1D/qvp_1D_'+date+'.nc', mode='w', format='NETCDF4')

############## After Praia - REDO August 29-31 #################
# Load data
file_dirs = glob.glob("/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4a/qvp_volume/2024083*")
file_dirs.append("/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4a/qvp_volume/20240829")

# loop through directories
for file_dir in file_dirs:
    # get the date from the directory name
    date = os.path.basename(file_dir)
    print(f"Processing data for {date}")

    file_paths = glob.glob(file_dir+"/*.nc")
    qvp = xr.open_mfdataset(file_paths, combine='by_coords')

    # Mask the data
    # missing data (blocking sectors)
    rainrate = qvp.RATE_CSU_BLENDED.where(qvp.RATE_CSU_BLENDED >=-30000, np.nan)
    dbz = qvp.DBZ.where(qvp.DBZ >=-30000, np.nan)

    #change -9999 missing data to zeros for rain, NaN for DBZ
    rainrate2 = rainrate.where(rainrate != -9999, 0)
    dbz2 = dbz.where(dbz != -9999, np.nan)

    # Average spatially at each level (results in time x height)
    rainrate_avg = rainrate2.mean(dim=['X', 'Y'],skipna=True)
    dbz_avg = dbz2.mean(dim=['X', 'Y'],skipna=True) 
    rainrate_cond = rainrate2.where(rainrate2>0).mean(dim=('X','Y'),skipna=True)

    # Save to netcdf file
    rainrate_avg.attrs['units'] = 'mm/h'
    rainrate_cond.attrs['units'] = 'mm/h'
    dbz_avg.attrs['units'] = 'dBZ'

    #combine into one dataset
    qvp1D= xr.Dataset({'rainrate_avg': rainrate_avg, 'dbz_avg': dbz_avg, 'rainrate_int': rainrate_cond})
    #qvp1D= xr.Dataset({'dbz_avg': dbz_avg})

    #add attributes
    qvp1D.attrs['title'] = 'Average rain rate and reflectivity at each altitude from SEA-POL QVP'
    qvp1D.attrs['description'] = 'Spatial means of rainrate from SEA-POL 45 degree elevation scans at different altitudes. Mean is the average over all valid data points (including zeros). Int (intensity) is the average over all data points where rainrate > 0.'
    qvp1D.attrs['source'] = 'SEA-POL Level4a QVP Volume'

    #save to netcdf
    qvp1D.to_netcdf('/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4a/qvp_1D/qvp_1D_'+date+'.nc', mode='w', format='NETCDF4')