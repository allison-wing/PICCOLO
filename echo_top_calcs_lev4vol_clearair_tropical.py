#%%
# This code is aimed for the VOL scans and doesn't include the gridded satellite.. 
##############  Import packages ############## 
import matplotlib.pyplot as plt
import matplotlib.colors as ticker
import matplotlib.colors as colors
import xarray as xr
import numpy as np
import os
import glob
import pandas as pd
import numpy as np
from datetime import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gs
from matplotlib.colors import ListedColormap
from datetime import datetime, timedelta



### Input the date that we are looping through
#start_time = datetime(2024, 8, 16, 8, 10, 0) 
#end_time   = datetime(2024, 9, 23, 16, 20, 59)  

start_time = datetime(2024, 9, 16, 0, 0, 0) 
end_time= datetime(2024, 9, 23, 16, 20, 59)

interval   = timedelta(minutes=10)
#%%

"""
Find the height (in km) of the highest X dBZ echo followed by a lower value or missing (-999/nan).
Parameters:
    reflectivity (xr.DataArray): 3D DataArray with dimensions ('z', 'y', 'x').
Returns:
    xr.DataArray: 2D DataArray with height in km of 20 dBZ echo, with lat/lon coords.
"""


# This function now accounts for clear air as well 
 
def find_height_DBZecho(DBZ_want, reflectivity, rainrate):
    """
    Finds the height index (Z) of the first occurrence of DBZ_want or higher
    from the top of the column. Works for data with (Z,Y,X) or (time,Z,Y,X),
    and always returns a DataArray with a time dimension.
    """

    # Ensure reflectivity has a time dimension
    if 'time' not in reflectivity.dims:
        reflectivity = reflectivity.expand_dims(time=[0])

    times = reflectivity.time

    all_heights = []
    all_confidence = []
    all_continuity = []
    all_gap_masks = []
    all_gap_thickness = []

    for t in times:
        # this is reflectivity
        ref_t = reflectivity.sel(time=t)
        z_dim, y_dim, x_dim = ref_t.shape

        heights = np.full((y_dim, x_dim), np.nan)
        confidence = np.full((y_dim, x_dim), np.nan)
        continuity = np.full((y_dim, x_dim), np.nan)

        # arrays for the additional clear air gap checks
        gap_mask = np.zeros((z_dim, y_dim, x_dim), dtype=bool)
        gap_thickness_arr = np.zeros((y_dim, x_dim))

        ref = ref_t.values  # for speed

        # this s for rain rate
        rr_t = rainrate.sel(time=t).isel(Z= 1) # selecting the time and z=1 km
        rr = rr_t.values # for speed

        for j in range(y_dim):
            for i in range(x_dim):
                column = ref[:, j, i]
                sfc_rr= rr[j, i]

                # 1. if there is a non-zero RR at surface 
                if sfc_rr > 0:
                     
                    
                    # loop from top down (just after 1 km ) ## previously was range(x_dim)
                    for z in reversed(range(1, z_dim)):
                        val = column[z] # value that at that level

                        if np.isnan(val):
                            continue
                        
                        # Check if reflectivity meets or exceeds threshold
                        if val >= DBZ_want:

                            heights[j, i] = z

                            # ---- CONFIDENCE MASK LOGIC ----
                            if z == z_dim - 1:
                                # no level above
                                confidence[j, i] = 0
                            else:
                                above = column[z + 1]

                                if np.isnan(above):
                                    confidence[j, i] = 0
                                elif np.isclose(above, -9999) or above < DBZ_want:
                                    confidence[j, i] = 1
                                else:  # above >= DBZ_want
                                    confidence[j, i] = 2
                            # --------------------------------
                            
                            # ---- ECHO CONTINUITY LOGIC ----
                            # Check if the echo is continous to the surface or there is another feature above 
                            gap_mask = np.zeros((z_dim, y_dim, x_dim), dtype=bool)

                            # select levels below the echo top (starting at 1km )
                            below = column[1:z]

                            # different conditions
                            is_nan = np.isnan(below)
                            is_clear = below == -9999  
                            is_real = (~is_nan) & (~is_clear)

                            # A. Fully continuous
                            if np.all(is_real):
                                continuity[j, i] = 1   # continuous to surface

                                gap_thickness_arr[j, i] = 0

                                

                            # B. Clear-air gap anywhere
                            elif np.any(is_clear):
                                continuity[j, i] = 0   # not continuous
                                
                                ### Additional checks 
                                # 1 identifying the actual gap levels
                                #gap_levels = np.where(is_clear)[0] + 1   # +1 because we sliced from 1
                                #gap_levels_arr[j, i] = gap_levels # storing

                                # absolute vertical indices
                                gap_indices = np.where(is_clear)[0] + 1

                                # store in 3D mask
                                gap_mask[gap_indices, j, i] = True

                                # store thickness
                                gap_thickness_arr[j, i] = len(gap_indices)

                            # C. Mix of real values and NaNs (but no clear air)
                            elif np.any(is_real) and np.any(is_nan):
                                continuity[j, i] = 2  # uncertain / missing data below

                                gap_thickness_arr[j, i] = 0

                            # --------------------------------

                            break

                # 2. NOT a valid rainrate... 
                else: 
                    is_nan = np.isnan(column)
                    is_clear = np.isclose(column, -9999)

                     # a. if all column is nan
                    if np.all(is_nan): 
                        heights[j,i] = np.nan
                        continue # because the array is already filled with nans originally 
                    
                    # b. if the all column is -9999, then CLEAR AIR
                    elif np.all(is_clear): 
                         heights[j, i] = -1
                         continue
                    
                    # c. if the column has an elevated dbz after 1 km, then VIRGA
                    elif np.any(column[2:] == DBZ_want):
                         heights[j, i] = -2
                         continue  
                    
                    # d.  mix of nans and -9999s
                    elif np.any(is_nan) and np.any(is_clear):
                        valid = is_nan | is_clear
                        # if majority is -9999, then lets set as clear air

                        if is_clear.sum() / valid.sum() >= 0.3:
                            heights[j, i] = -1

                        # if not just leave as nan
                        else:
                            heights[j, i] = np.nan
 
                   
        # store results for this time step
        height_da = xr.DataArray(
            heights,
            dims=("Y", "X"),
            coords={"Y": ref_t.Y, "X": ref_t.X},
        )
        all_heights.append(height_da)

        conf_da = xr.DataArray(
            confidence,
            dims=("Y", "X"),
            coords={"Y": ref_t.Y, "X": ref_t.X},
            name=f"echo_top_{DBZ_want}dBZ_confidence",
        )
        all_confidence.append(conf_da)

        continuity_da = xr.DataArray(
            continuity,
            dims=("Y", "X"),
            coords={"Y": ref_t.Y, "X": ref_t.X},
            name=f"echo_top_{DBZ_want}dBZ_continuity",
        )
        all_continuity.append(continuity_da)

        gap_mask_da = xr.DataArray(
            gap_mask,
            dims=("Z", "Y", "X"),
            coords={"Z": ref_t.Z, "Y": ref_t.Y, "X": ref_t.X},
            name=f"echo_top_{DBZ_want}dBZ_gap_mask",
        )
        all_gap_masks.append(gap_mask_da)

        gap_thickness_da = xr.DataArray(
            gap_thickness_arr,
            dims=("Y", "X"),
            coords={"Y": ref_t.Y, "X": ref_t.X},
            name=f"echo_top_{DBZ_want}dBZ_gap_thickness",
        )
        all_gap_thickness.append(gap_thickness_da)


    # combine across time dimension
    height_all = xr.concat(all_heights, dim=times).assign_coords(time=times)
    height_all.name = f"echo_top_{DBZ_want}dBZ"
    height_all.attrs = {
        "description": f"Height index (Z) of {DBZ_want} dBZ echo top, clear air= -1, virga= -2",
        "values": "-1, -2, 0+",
        "units": "index of vertical level",

    }
    conf_all = xr.concat(all_confidence, dim=times).assign_coords(time=times)
    conf_all.name = f"echo_top_{DBZ_want}dBZ_confidence"
    conf_all.attrs = {
        "description": (
            "Confidence mask based on reflectivity one level above identified echo top: "
            "0=NaN above (maybe), 1=clear or < threshold (certain), 2=reflective above (special)"
        ),
        "values": "0,1,2",
    }

    continuity_all = xr.concat(all_continuity, dim=times).assign_coords(time=times)
    continuity_all.name = f"echo_top_{DBZ_want}dBZ_continuity"
    continuity_all.attrs = {
        "description": (
            "Continuity mask indicating if echo is continuous to surface: "
            "0=not continuous (clear air gap), 1=continuous, 2=uncertain/missing data below"
        ),
        "values": "0,1,2",
    }

    gap_mask_all = xr.concat(all_gap_masks, dim=times).assign_coords(time=times)
    gap_thickness_all = xr.concat(all_gap_thickness, dim=times).assign_coords(time=times)

    return height_all, conf_all, continuity_all, gap_mask_all, gap_thickness_all


# max reflectivity  

def find_max_DBZ(reflectivity):
    """
    Finds the maximum reflectivity (dBZ) in each vertical column.
    Works with data shaped as (Z, Y, X) or (time, Z, Y, X).
    Always returns a DataArray with a time dimension.
    """

    # --- Case 1: If there's no time dimension, add a dummy one ---
    if 'time' not in reflectivity.dims:
        reflectivity = reflectivity.expand_dims(time=[0])

    # --- Compute maximum along the vertical (Z) dimension ---
    # (skip NaNs automatically)
    max_dbz = reflectivity.max(dim='Z', skipna=True)

    # --- Ensure attributes and name ---
    max_dbz.name = 'max_reflectivity'
    max_dbz.attrs.update({
        'description': 'Maximum reflectivity (dBZ) in each vertical column',
        'long_name': 'Column maximum reflectivity',
        'units': 'dBZ'
    })

    return max_dbz

def keeping_120km(echo_var): 
    # keeping just within 120km range 
        
    lat0 = echo_var.latitude.sel(X=0, Y=0, method="nearest").item()
    lon0 = echo_var.longitude.sel(X=0, Y=0, method="nearest").item()

    # calculate distance 
    dlat_km = (echo_var.latitude - lat0) * 111.0
    dlon_km = (echo_var.longitude - lon0) * 111.0 * np.cos(np.deg2rad(lat0))

    distance_km = np.sqrt(dlat_km**2 + dlon_km**2)

    # subset data based on distance 
    height_120km = echo_var.where(distance_km <= 120)

    return height_120km

# not found count
nf_count= 1
#%%
# Generate all time points (every 10 minutes)
time_steps = []
current_time = start_time
while current_time <= end_time:
    time_steps.append(current_time)
    current_time += interval

datasets= []

#### Open the level 4 VOL file
path= '/bell-scratch2/seapol/DATA/PICCOLO/qc_data/level4/uncompressed/'
file= 'PICCOLO_level4_volume_3D.nc'
rad_path= path+ file
rad_all= xr.open_dataset(rad_path)

time_index = rad_all.get_index("time")
rad_all = rad_all.isel(time=~rad_all.get_index("time").duplicated()) # remove duplicated times if any

times_xr = rad_all['time'].values # get the file times 

####################################
for dt in time_steps:

        formatted = dt.strftime("%Y%m%d_%H%M")
        date= dt.strftime("%Y%m%d")
        time_dt = pd.to_datetime(formatted, format='%Y%m%d_%H%M')
        time_np = np.datetime64(time_dt, 'ns')
        
        try:
            # if we don't have that time in dataset... 
            if not np.isin(time_np, times_xr):
                raise ValueError(f"Time {formatted} ({time_np}) not found in dataset.")

        except (ValueError):
                
                nf_count= nf_count + 1
                continue  # Skip this iteration

        #### Select the running time ...
        rad = rad_all.sel(time=[time_np]) 

        # file found - proceed to calculations
        ##########################
        merged= rad

        # Get numpy array of datetime64[ns]
        dt_array = merged.start_time.values  # e.g. numpy.datetime64[ns]
        dt_minutes = dt_array.astype('datetime64[m]') # 1. Convert to minutes since epoch
        minutes_int = dt_minutes.astype('int')# 2. Get integer minutes since epoch
        rounded_minutes = ((minutes_int + 5) // 10) * 10  # add 5 for rounding to nearest # 3. Round to nearest 10 minutes
        rounded_dt = rounded_minutes.astype('datetime64[m]') # 4. Convert back to datetime64[m]
        merged['time'] = rounded_dt.astype('datetime64[ns]') # 5. Store it (you can cast to ns if needed)

        # Calculating echo top height and assigning to array.... 
        reflectivity= merged.DBZ
        rainrate= merged.RAINRATE

        merged['RAINRATE_zero'] = merged['RAINRATE'].where(merged['RAINRATE'] != -9999, other=0) # 'RAINRATE' with level 4b # RATE_CSU_BLENDED
        #merged['echo_height_neg5']= find_height_DBZecho(-5, reflectivity)
        #merged['echo_height_0']= find_height_DBZecho(0, reflectivity, rainrate)
        #merged['echo_height_5']= find_height_DBZecho(5, reflectivity, rainrate)
        merged['echo_height_10'], merged['echo_10_conf'], merged['echo_10_continuity'], merged['clear_gap_mask'], merged['clear_gap_thickness'] = find_height_DBZecho(10, reflectivity, rainrate)
        #merged['echo_height_15']= find_height_DBZecho(15, reflectivity, rainrate)
        #merged['echo_height_20']= find_height_DBZecho(20, reflectivity, rainrate)
        #merged['echo_height_30']= find_height_DBZecho(30, reflectivity, rainrate)
        #merged['echo_height_40']= find_height_DBZecho(40, reflectivity, rainrate) # Not running the other XdBZs rihg now, just 10 dBZ
        merged['max_dBZ']= find_max_DBZ(reflectivity)

        # Keep only data within 120 km
        merged['echo_height_10']= keeping_120km(merged['echo_height_10'])
        merged['echo_10_conf']= keeping_120km(merged['echo_10_conf'])
        merged['echo_10_continuity']= keeping_120km(merged['echo_10_continuity'])
        merged['clear_gap_mask']= keeping_120km(merged['clear_gap_mask'])
        merged['clear_gap_thickness']= keeping_120km(merged['clear_gap_thickness'])

        # Append each xarray... 
        datasets.append(merged)

        # Convert X and Y from meters to kilometers
        merged_km = merged.assign_coords(
            X = merged.X / 1000.0,
            Y = merged.Y / 1000.0
        )

        merged_km.X.attrs['units'] = 'km'
        merged_km.Y.attrs['units'] = 'km'

        # export the individual files...
        out_file= 'radsat_calc_'+ formatted + '.nc'
        out_path = '/bell-scratch/deliancb/piccolo_mpi_proj/output_nc/echo_top_calcs_VOL_feb/' + out_file

        # Check if file exists and delete it 
        if os.path.exists(out_path):
            print(f"{out_path} already exists. Removing it...")
            os.remove(out_path)

        # Now write the NetCDF
        merged_km.to_netcdf(out_path)
        print(f"Saved {out_path}")


# If you want a concat file then uncomment below line
#ds_combined = xr.concat(datasets, dim="time")
print ('Number of files not found (times where there was no match): '+ str(nf_count))


