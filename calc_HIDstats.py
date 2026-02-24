## Calculate echo base and top height for each frame of SEA-POL volume data 

import xarray as xr
import numpy as np
import pandas as pd

# read in data
seapol = xr.open_dataset('/huracan/tank4/cornell/ORCESTRA/sea-pol/qc_data/level4v1.2/PICCOLO_level4_volume_3D.nc')


# Make regular 10-minute time series
#start_time = np.datetime64('2024-08-16T08:00:00')
#start_time = np.datetime64('2024-09-13T00:00:00') #East/west division
start_time = np.datetime64('2024-08-28T20:20:00') #AP start
end_time = np.datetime64('2024-09-23T16:50:00')
time10m = pd.date_range(start_time, end_time, freq='10 min')
time10m = pd.to_datetime(time10m)

# Define HID bins for histogram
hid_bins = np.arange(0.5,11.5,1)

# Pre-define arrays to store results, len(time10m) rows, len(hid_bins) colums
all_hid_hist = np.full((len(time10m),len(hid_bins)-1), np.nan)  # -1 because histogram returns one fewer bin edges than bins
all_hid_density = np.full((len(time10m),len(hid_bins)-1), np.nan)  # -1 because histogram returns one fewer bin edges than bins

# Loop over times
#for i in range(0, np.size(seapol.time[:-1])): #loop over all times
for i in range(0,len(time10m)):  # loop over all times in the 10-minute series (VOL1 only)

    #check if there is data for this time
    #if seapol.time[i].values not in seapol.time.values:
    if time10m[i] not in seapol.time.values: #VOL1 only
        #print(f"No data for time {seapol.time[i].values}")
        print(f"No data for time {time10m[i]}") #VOL1 only
        hid_hist = np.full(len(hid_bins)-1, np.nan)  # -1 because histogram returns one fewer bin edges than bins
        hid_density = np.full(len(hid_bins)-1, np.nan)  # -1 because histogram returns one fewer bin edges than bins
    else:
        #print('Time:', seapol.time[indexAP + i].values)
        #map_dbz = seapol.DBZ[indexAP+i,:,:,:]
         
        print('Time:', seapol.time.sel(time=time10m[i]).values) #VOL1 only
        vol_HID = seapol.HID_CSU.sel(time=time10m[i]) #VOL1 only
        vol_HID = vol_HID.where(vol_HID !=-9999, np.nan) #set clear air echoes to NaN (no data possible already NaN)

        # Calculate counts and density for this scene
        hid_hist,bin_edges = np.histogram(vol_HID.values.flatten(),bins=hid_bins,density=False)
        hid_density,bin_edges = np.histogram(vol_HID.values.flatten(),bins=hid_bins, density=True)
        
    # outside if statement so that if there is no data for this time, the histogram and density will be NaN and will be stored as such in the global arrays
    all_hid_hist[i,:] = hid_hist
    all_hid_density[i,:] = hid_density
    

# Save the dataset to netcdf
ds = xr.Dataset(data_vars={'hid_hist': (['time', 'hid_bins'], all_hid_hist),
                'hid_density': (['time', 'hid_bins'], all_hid_density),'bin_edges':(['hid_bin_edges'],bin_edges)},
                coords={'time': time10m,
                        'hid_bins': bin_edges[:-1],
                        'hid_bin_edges':bin_edges})
ds.to_netcdf('../../data/SEA-POLv1.2_vol_HID_hist_density.nc')