# %% [markdown]
# # Spatial plots of CWV and CAPE and shear

# %% [markdown]
# From ERA-5

# %%
import geopy.distance
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import gridspec
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import seaborn as sns
import numpy as np
import xarray as xr
import json
import cartopy
import pandas as pd
from datetime import datetime, timedelta
import cftime
#from adjustText import adjust_text
import glob
import os

# %% [markdown]
# ### Set region to plot and start and end time

# %%
lonMin, lonMax = -62.1, -9.9
latMin, latMax = -2.1, 22.1
#latMin, latMax = -10.1,30.1

start_mon = 'Aug.'
start_day = '10'
start_time = '08-'+start_day+'T00:00:00'
end_mon = 'Sep.'
end_day = '30'
end_time = '09-'+end_day+'T00:00:00'

# %% [markdown]
# ### Plot CWV

# %%
# For each year, spatial plot of CWV

fig = plt.figure(figsize = (20, 12))
gs = gridspec.GridSpec(6, 5, width_ratios=[1, 1, 1, 1, 1])

# Filename base for ERA-5 CWV data
filebase_CWV = "/huracan/tank4/cornell/ORCESTRA/era5/total_column_water_vapour/"

years = np.arange(1996,2025)  # 1998-2024
iyear = 0
for yy in years:
    # print year
    print("Processing year:", yy)
    
    # Load ERA-5 CWV data
    file_paths_CWV_08 = glob.glob(filebase_CWV + str(yy) + '08/' + "*.nc")
    file_paths_CWV_09 = glob.glob(filebase_CWV + str(yy) + '09/' + "*.nc")
    file_paths_CWV = file_paths_CWV_08 + file_paths_CWV_09
    era5_cwv = xr.open_mfdataset(file_paths_CWV,combine='by_coords')
    
    #Extract for time period of campaign and region
    cwv_yy = era5_cwv.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin,lonMax),valid_time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))

    #Take mean over campaign
    cwv_yy_mean = cwv_yy.mean(dim='valid_time')
    
    # Create plot
    ax1 = fig.add_subplot(gs[iyear],projection=ccrs.PlateCarree())
    ax1.coastlines(resolution = '50m',alpha=0.5)
    #ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
    gl.top_labels = False
    gl.left_labels = False
    if iyear==25 or iyear==26 or iyear==27 or iyear==28:
        print(iyear)    
        gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.xlabel_style = {'size': 14}
    else: 
        gl.bottom_labels = False
    if iyear==4 or iyear==9 or iyear==14 or iyear==19 or iyear==24 or iyear==28:
        print(iyear)
        gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
        gl.yformatter = LATITUDE_FORMATTER
        gl.ylabel_style = {'size': 14}
    else:
        gl.right_labels = False
    
    # Plot this years CWV as filled contours
    plt.contourf(cwv_yy_mean.longitude.values, 
                 cwv_yy_mean.latitude.values,
                 cwv_yy_mean.tcwv.values,
                 cmap = 'viridis',levels = np.arange(26, 60, 4),extend='both')

    cbar=plt.colorbar(location='left',orientation='vertical',shrink = 0.75, pad=0.02)
    cbar.set_label('CWV (kg m$^{-2}$)',fontsize=14)
    cbar.set_ticks(np.arange(26, 60, 4))
    cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(26, 60, 4)],fontsize=14)

    # Set the plot limits and title
    ax1.set_title(str(yy),size=14)
    ax1.set_xlim(lonMin, lonMax)
    ax1.set_ylim(latMin, latMax)
    
    era5_cwv.close()
    cwv_yy.close()
    iyear=iyear+1
    
sns.set_context('paper') 
plt.tight_layout()  # Adjust subplots to fit into figure area.
plt.savefig('Figure_ERA5_CWV_spatial.png', dpi=300)

# %%
# For each year, spatial plot of CAPE

fig = plt.figure(figsize = (20, 12))
gs = gridspec.GridSpec(6, 5, width_ratios=[1, 1, 1, 1, 1])

# Filename base for ERA-5 CAPE data
filebase_CAPE = "/huracan/tank4/cornell/ORCESTRA/era5/convective_available_potential_energy/"

years = np.arange(1996,2025)  # 1998-2024
iyear = 0
for yy in years:
    # print year
    print("Processing year:", yy)
    
    # Load ERA-5 CAPE data
    file_paths_CAPE_08 = glob.glob(filebase_CAPE + str(yy) + '08/' + "*.nc")
    file_paths_CAPE_09 = glob.glob(filebase_CAPE + str(yy) + '09/' + "*.nc")
    file_paths_CAPE = file_paths_CAPE_08 + file_paths_CAPE_09
    era5_cape = xr.open_mfdataset(file_paths_CAPE,combine='by_coords')

    #Extract for time period of campaign and region
    cape_yy = era5_cape.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin,lonMax),valid_time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))

    #Take mean over campaign
    cape_yy_mean = cape_yy.mean(dim='valid_time')

    # Create plot
    ax1 = fig.add_subplot(gs[iyear],projection=ccrs.PlateCarree())
    ax1.coastlines(resolution = '50m',alpha=0.5)
    #ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
    gl.top_labels = False
    gl.left_labels = False
    if iyear==25 or iyear==26 or iyear==27 or iyear==28:
        print(iyear)    
        gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.xlabel_style = {'size': 14}
    else: 
        gl.bottom_labels = False
    if iyear==4 or iyear==9 or iyear==14 or iyear==19 or iyear==24 or iyear==28:
        print(iyear)
        gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
        gl.yformatter = LATITUDE_FORMATTER
        gl.ylabel_style = {'size': 14}
    else:
        gl.right_labels = False

    # Plot this years CAPE as filled contours
    plt.contourf(cape_yy_mean.longitude.values, 
                 cape_yy_mean.latitude.values,
                 cape_yy_mean.cape.values,
                 cmap = 'viridis',levels=np.arange(200,2000,300),extend='both')

    cbar=plt.colorbar(location='left',orientation='vertical',shrink = 0.75, pad=0.02)
    cbar.set_label('CAPE (K kg$^{-1}$)',fontsize=14)
    cbar.set_ticks(np.arange(200, 2000, 300))
    cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(200, 2000, 300)],fontsize=14)

    # Set the plot limits and title
    ax1.set_title(str(yy),size=14)
    ax1.set_xlim(lonMin, lonMax)
    ax1.set_ylim(latMin, latMax)

    era5_cape.close()
    cape_yy.close()
    iyear=iyear+1
    
sns.set_context('paper') 
plt.tight_layout()  # Adjust subplots to fit into figure area.
plt.savefig('Figure_ERA5_CAPE_spatial.png', dpi=300)

# %% [markdown]
# ### Plot shear

# %% [markdown]
# First plot just 2024 shear

# %%
# 10-m winds
# Filename base
filebase = "/huracan/tank4/cornell/ORCESTRA/era5/10m_u_component_of_wind/2024"

file_paths08 = glob.glob(filebase + '08/' + "*.nc")
file_paths09 = glob.glob(filebase + '09/' + "*.nc")
file_paths = file_paths08 + file_paths09

# open the data
era5_usfc= xr.open_mfdataset(file_paths, combine='by_coords')

filebase = "/huracan/tank4/cornell/ORCESTRA/era5/10m_v_component_of_wind/2024"

file_paths08 = glob.glob(filebase + '08/' + "*.nc")
file_paths09 = glob.glob(filebase + '09/' + "*.nc")
file_paths = file_paths08 + file_paths09

# open the data
era5_vsfc = xr.open_mfdataset(file_paths, combine='by_coords')

# %%
# winds on pressure levels
# Filename base
filebase = "/huracan/tank4/cornell/ORCESTRA/era5/u_component_of_wind/2024"

file_paths08 = glob.glob(filebase + '08/' + "*.nc")
file_paths09 = glob.glob(filebase + '09/' + "*.nc")
file_paths = file_paths08 + file_paths09

# open the data
era5_u= xr.open_mfdataset(file_paths, combine='by_coords')

# Filename base
filebase = "/huracan/tank4/cornell/ORCESTRA/era5/v_component_of_wind/2024"

file_paths08 = glob.glob(filebase + '08/' + "*.nc")
file_paths09 = glob.glob(filebase + '09/' + "*.nc")
file_paths = file_paths08 + file_paths09

# open the data
era5_v= xr.open_mfdataset(file_paths, combine='by_coords')


# %%
# calculate wind shear

# surface to 600 hPa [approx equivalent of 0-6 km]
era5_shear06 = np.sqrt( (era5_u.u.sel(pressure_level=600) - era5_usfc.u10)**2 + (era5_v.v.sel(pressure_level=600) - era5_vsfc.v10)**2 )

# surface to 850 hPa [approx equivalent of 0-2 km]
era5_shear02 = np.sqrt( (era5_u.u.sel(pressure_level=850) - era5_usfc.u10)**2 + (era5_v.v.sel(pressure_level=850) - era5_vsfc.v10)**2 )

# %%
# Take mean over campaign time period
shear06 = era5_shear06.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin,lonMax),valid_time=slice('2024-'+start_time,'2024-'+end_time)).mean(dim='valid_time')
shear02 = era5_shear02.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin,lonMax),valid_time=slice('2024-'+start_time,'2024-'+end_time)).mean(dim='valid_time')


# %%
plt.figure(figsize = (18, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution = '50m',alpha=0.5)
#ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 18}
gl.ylabel_style = {'size': 18}

plt.contourf(shear06.longitude.values, shear06.latitude.values,shear06.values, cmap = 'viridis', levels = np.arange(0,22,2), extend='both')

# Set the plot limits and title
plt.title(start_day+' '+start_mon+' - '+end_day+' '+end_mon+' 2024 Mean Sfc - 600 hPa Shear (m/s)',size=18)
ax.set_xlim(lonMin, lonMax)
ax.set_ylim(latMin, latMax) 

# Add colorbar
cbar=plt.colorbar(location='right',orientation='vertical',shrink = 0.9, pad=0.02)
cbar.set_label('S06 (m s$^{-1}$)',fontsize=14)
cbar.set_ticks(np.arange(0, 22, 2))
cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(0, 22, 2)],fontsize=14)  

plt.savefig('Figure_ERA5_Shear06_spatial_2024.png', dpi=300)

# %%
plt.figure(figsize = (18, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution = '50m',alpha=0.5)
#ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 18}
gl.ylabel_style = {'size': 18}

plt.contourf(shear02.longitude.values, shear02.latitude.values,shear02.values, cmap = 'viridis', levels = np.arange(0,12,2), extend='both')

# Set the plot limits and title
plt.title(start_day+' '+start_mon+' - '+end_day+' '+end_mon+' 2024 Mean Sfc - 850 hPa Shear (m/s)',size=18)
ax.set_xlim(lonMin, lonMax)
ax.set_ylim(latMin, latMax) 

# Add colorbar
cbar=plt.colorbar(location='right',orientation='vertical',shrink = 0.9, pad=0.02)
cbar.set_label('S02 (m s$^{-1}$)',fontsize=14)
cbar.set_ticks(np.arange(0, 12, 2))
cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(0, 12, 2)],fontsize=14)  

plt.savefig('Figure_ERA5_Shear02_spatial_2024.png', dpi=300)

# %%
# For each year, spatial plot of sfc to 600 hPa shear

fig = plt.figure(figsize = (20, 12))
gs = gridspec.GridSpec(6, 5, width_ratios=[1, 1, 1, 1, 1])

# Filename base for ERA-5 10m wind data
filebase_10m = "/mars/tank3/era5/uv_10m_wind/"

# Load ERA-5 u data because all years are in one file
era5_u600 = xr.open_dataset('/mars/tank3/era5/pressure_variables/u_component_of_wind/u_600.nc')

# Load ERA-5 v data because all years are in one file
era5_v600 = xr.open_dataset('/mars/tank3/era5/pressure_variables/v_component_of_wind/v_600.nc')

years = np.arange(1990,2019)
iyear = 0
for yy in years:
    # print year
    print("Processing year:", yy)
    
    # Load ERA-5 10m wind data
    file_paths_10m_08 = glob.glob(filebase_10m + str(yy) + '-08.nc')
    file_paths_10m_09 = glob.glob(filebase_10m + str(yy) + '-09.nc')
    file_paths_10m = file_paths_10m_08 + file_paths_10m_09
    era5_10m = xr.open_mfdataset(file_paths_10m,combine='by_coords')

    #Extract for time period of campaign and region
    era5_10m_yy = era5_10m.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin+360,lonMax+360),time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))
    u600_yy = era5_u600.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin+360,lonMax+360),time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))
    v600_yy = era5_v600.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin+360,lonMax+360),time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))

    # Calculate wind shear
    # surface to 600 hPa (approx equivalent of 0-6 km)
    shear06_yy = np.sqrt((u600_yy.u - era5_10m_yy.u10)**2 + (v600_yy.v - era5_10m_yy.v10)**2)

    #Take mean over campaign
    shear06_mean = shear06_yy.mean(dim='time')

    # Create plot
    ax1 = fig.add_subplot(gs[iyear],projection=ccrs.PlateCarree())
    ax1.coastlines(resolution = '50m',alpha=0.5)
    #ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
    gl.top_labels = False
    gl.left_labels = False
    if iyear==25 or iyear==26 or iyear==27 or iyear==28 or iyear==29:
        print(iyear)    
        gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.xlabel_style = {'size': 14}
    else: 
        gl.bottom_labels = False
    if iyear==4 or iyear==9 or iyear==14 or iyear==19 or iyear==24 or iyear==29:
        print(iyear)
        gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
        gl.yformatter = LATITUDE_FORMATTER
        gl.ylabel_style = {'size': 14}
    else:
        gl.right_labels = False

    # Plot this years shear as filled contours
    plt.contourf(shear06_mean.longitude.values, 
                 shear06_mean.latitude.values,
                 shear06_mean.values,
                 cmap = 'viridis',levels=np.arange(0,22,2),extend='both')

    cbar=plt.colorbar(location='left',orientation='vertical',shrink = 0.75, pad=0.02)
    cbar.set_label('S06 (m s$^{-1}$)',fontsize=14)
    cbar.set_ticks(np.arange(0, 22, 2))
    cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(0, 22, 2)],fontsize=14)

    # Set the plot limits and title
    ax1.set_title(str(yy),size=14)
    ax1.set_xlim(lonMin, lonMax)
    ax1.set_ylim(latMin, latMax)

    era5_10m.close()
    era5_10m_yy.close()
    u600_yy.close()
    v600_yy.close()
    u850_yy.close()
    v850_yy.close()
    iyear=iyear+1

# Add 2024
yy = 2024    
# Create plot
ax1 = fig.add_subplot(gs[iyear],projection=ccrs.PlateCarree())
ax1.coastlines(resolution = '50m',alpha=0.5)
#ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
gl = ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
gl.top_labels = False
gl.left_labels = False
if iyear==25 or iyear==26 or iyear==27 or iyear==29:
    print(iyear)    
    gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 14}
else: 
    gl.bottom_labels = False
if iyear==4 or iyear==9 or iyear==14 or iyear==19 or iyear==24 or iyear==29:
    print(iyear)
    gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'size': 14}
else:
    gl.right_labels = False

# Plot this years shear as filled contours
plt.contourf(shear06.longitude.values, 
            shear06.latitude.values,
            shear06.values,
            cmap = 'viridis',levels=np.arange(0,22,2),extend='both')

cbar=plt.colorbar(location='left',orientation='vertical',shrink = 0.75, pad=0.02)
cbar.set_label('S06 (m s$^{-1}$)',fontsize=14)
cbar.set_ticks(np.arange(0, 22, 5))
cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(0, 22, 5)],fontsize=14)

# Set the plot limits and title
ax1.set_title(str(yy),size=14)
ax1.set_xlim(lonMin, lonMax)
ax1.set_ylim(latMin, latMax)    

sns.set_context('paper') 
plt.tight_layout()  # Adjust subplots to fit into figure area.
plt.savefig('Figure_ERA5_Shear06_spatial.png', dpi=300)

# %%
# For each year, spatial plot of sfc to 850 hPa shear

fig = plt.figure(figsize = (20, 12))
gs = gridspec.GridSpec(6, 5, width_ratios=[1, 1, 1, 1, 1])

# Filename base for ERA-5 10m wind data
filebase_10m = "/mars/tank3/era5/uv_10m_wind/"

# Load ERA-5 u data because all years are in one file
era5_u850 = xr.open_dataset('/mars/tank3/era5/pressure_variables/u_component_of_wind/u_850.nc')

# Load ERA-5 v data because all years are in one file
era5_v850 = xr.open_dataset('/mars/tank3/era5/pressure_variables/v_component_of_wind/v_850.nc')

years = np.arange(1990,2019)
iyear = 0
for yy in years:
    # print year
    print("Processing year:", yy)
    
    # Load ERA-5 10m wind data
    file_paths_10m_08 = glob.glob(filebase_10m + str(yy) + '-08.nc')
    file_paths_10m_09 = glob.glob(filebase_10m + str(yy) + '-09.nc')
    file_paths_10m = file_paths_10m_08 + file_paths_10m_09
    era5_10m = xr.open_mfdataset(file_paths_10m,combine='by_coords')

    #Extract for time period of campaign and region
    era5_10m_yy = era5_10m.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin+360,lonMax+360),time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))
    u850_yy = era5_u850.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin+360,lonMax+360),time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))
    v850_yy = era5_v850.sel(latitude=slice(latMax,latMin),longitude=slice(lonMin+360,lonMax+360),time=slice(str(yy)+'-'+start_time,str(yy)+'-'+end_time))

    # Calculate wind shear
    # surface to 850 hPa (approx equivalent of 0-2 km)
    shear02_yy = np.sqrt((u850_yy.u - era5_10m_yy.u10)**2 + (v850_yy.v - era5_10m_yy.v10)**2)

    #Take mean over campaign
    shear02_mean = shear02_yy.mean(dim='time')

    # Create plot
    ax1 = fig.add_subplot(gs[iyear],projection=ccrs.PlateCarree())
    ax1.coastlines(resolution = '50m',alpha=0.5)
    #ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
    gl.top_labels = False
    gl.left_labels = False
    if iyear==25 or iyear==26 or iyear==27 or iyear==28 or iyear==29:
        print(iyear)    
        gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.xlabel_style = {'size': 14}
    else: 
        gl.bottom_labels = False
    if iyear==4 or iyear==9 or iyear==14 or iyear==19 or iyear==24 or iyear==29:
        print(iyear)
        gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
        gl.yformatter = LATITUDE_FORMATTER
        gl.ylabel_style = {'size': 14}
    else:
        gl.right_labels = False

    # Plot this years shear as filled contours
    plt.contourf(shear02_mean.longitude.values, 
                 shear02_mean.latitude.values,
                 shear02_mean.values,
                 cmap = 'viridis',levels=np.arange(0,12,2),extend='both')

    cbar=plt.colorbar(location='left',orientation='vertical',shrink = 0.75, pad=0.02)
    cbar.set_label('S02 (m s$^{-1}$)',fontsize=14)
    cbar.set_ticks(np.arange(0, 12, 2))
    cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(0, 12, 2)],fontsize=14)

    # Set the plot limits and title
    ax1.set_title(str(yy),size=14)
    ax1.set_xlim(lonMin, lonMax)
    ax1.set_ylim(latMin, latMax)

    era5_10m.close()
    era5_10m_yy.close()
    u600_yy.close()
    v600_yy.close()
    u850_yy.close()
    v850_yy.close()
    iyear=iyear+1

# Add 2024
yy = 2024    
# Create plot
ax1 = fig.add_subplot(gs[iyear],projection=ccrs.PlateCarree())
ax1.coastlines(resolution = '50m',alpha=0.5)
#ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha = 0.25)
gl = ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,alpha = 0.25)
gl.top_labels = False
gl.left_labels = False
if iyear==25 or iyear==26 or iyear==27 or iyear==29:
    print(iyear)    
    gl.xlocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 14}
else: 
    gl.bottom_labels = False
if iyear==4 or iyear==9 or iyear==14 or iyear==19 or iyear==24 or iyear==29:
    print(iyear)
    gl.ylocator = mticker.FixedLocator([0, 5, 10, 15, 20])
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'size': 14}
else:
    gl.right_labels = False

# Plot this years shear as filled contours
plt.contourf(shear02.longitude.values, 
            shear02.latitude.values,
            shear02.values,
            cmap = 'viridis',levels=np.arange(0,12,2),extend='both')

cbar=plt.colorbar(location='left',orientation='vertical',shrink = 0.75, pad=0.02)
cbar.set_label('S02 (m s$^{-1}$)',fontsize=14)
cbar.set_ticks(np.arange(0, 12, 2))
cbar.set_ticklabels([f'{x:.0f}' for x in np.arange(0, 12, 2)],fontsize=14)

# Set the plot limits and title
ax1.set_title(str(yy),size=14)
ax1.set_xlim(lonMin, lonMax)
ax1.set_ylim(latMin, latMax)    

sns.set_context('paper') 
plt.tight_layout()  # Adjust subplots to fit into figure area.
plt.savefig('Figure_ERA5_Shear02_spatial.png', dpi=300)

