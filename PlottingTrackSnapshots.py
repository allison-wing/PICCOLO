#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#Open the files
regdata = xr.open_dataset('ERA-5_CampaignTrack_Regular_Variables_2016.nc')
budgdata = xr.open_dataset('ERA-5_CampaignTrack_Budget_Variables_2016.nc')

#Plot a variable of your choice
#Chose precip and converting from m/s to mm/hr
precipdata = regdata.TotPrecipratebox * 1000 * 3600
#Set the bound of colorbar based on max and min values of data
vmax = np.nanmax(precipdata)
vmin = np.nanmin(precipdata)
#Loop through all the track points
for t in range(0,len(precipdata)):
    #Set up the plot for the snapshot (vmin and vmax are colorbar ranges, extent is setting the numbers of x and y ticks)
    plt.imshow(precipdata[t], vmin=vmin, vmax=vmax, extent=[-1.25,1.25,-1.25,1.25], cmap='YlGn')
    plt.xlabel("Longitude [Degrees]", fontweight='bold')
    plt.ylabel("Latitude [Degrees]", fontweight='bold')
    plt.title("Precipitation Rate Following \nWithin 1.25 Degrees of Track", fontweight='bold')
    cbar = plt.colorbar()
    cbar.set_label("Precip. Rate [mm/hr]", fontweight='bold')
    cbar.ax.set_yticklabels(cbar.ax.get_yticks(),weight='bold')
    plt.gca().set_xticklabels(plt.gca().get_xticks(),weight='bold')
    plt.gca().set_yticklabels(plt.gca().get_yticks(),weight='bold')
    plt.show()
    txt=str(t+1)
    plt.savefig('/home/awing/orcestra/precipratesnaps/FieldCampaign_PrecipRate_Snapshot_'+txt.zfill(3)+'.png')
    plt.close()
