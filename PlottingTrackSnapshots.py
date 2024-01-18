#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#year of interest
year = 2016

#Open the files
regdata = xr.open_dataset('ERA-5_CampaignTrack_Regular_Variables_'+str(year)+'.nc')
budgdata = xr.open_dataset('ERA-5_CampaignTrack_Budget_Variables_'+str(year)+'.nc')

#Grab time/track info
mon = np.array(regdata.month)
mon = mon.astype("int")
day = np.array(regdata.day)
day = day.astype("int")
hour = np.array(regdata.hour)
hour = hour.astype("int")

lat = regdata.latitude
lon = regdata.longitude-360
#Chose data field to plot
data = regdata.TotPrecipratebox * 1000 * 3600 #convert from m/s to mm/hr
#data = -1*regdata.OLR

#Set the bound of colorbar based on max and min values of data
vmax = np.nanmax(data)
vmin = np.nanmin(data)

#Loop through all the track points
for t in range(0,len(data)):
    #Set up the plot for the snapshot (vmin and vmax are colorbar ranges, extent is setting the numbers of x and y ticks)
    #plt.imshow(data[t], vmin=vmin, vmax=vmax, extent=[-1.25,1.25,-1.25,1.25], cmap='gray') #YlGn #gray
    plt.imshow(data[t], vmin=vmin, vmax=vmax, extent=[np.nanmin(lon[t]),np.nanmax(lon[t]),np.nanmin(lat[t]),np.nanmax(lat[t])], cmap='YlGn')
    plt.xlabel("Longitude [Degrees]", fontweight='bold')
    plt.ylabel("Latitude [Degrees]", fontweight='bold')
    #plt.title("Precip. \nWithin 1.25 Degrees of Track", fontweight='bold')
    #plt.title('OLR: ' + str(year)+'-'+"{0:0=2d}".format(mon[t])+'-'+"{0:0=2d}".format(day[t])+'-'+str(hour[t])+'Z \n centered around Meteor track', fontweight='bold') 
    plt.title('Precip: ' + str(year)+'-'+"{0:0=2d}".format(mon[t])+'-'+"{0:0=2d}".format(day[t])+'-'+str(hour[t])+'Z \n centered around Meteor track', fontweight='bold')
    cbar = plt.colorbar()
    cbar.set_label("Precip. Rate [mm/hr]", fontweight='bold')
    #cbar.set_label("OLR [W//m^2]", fontweight='bold')
    cbar.ax.set_yticklabels(cbar.ax.get_yticks(),weight='bold')
    plt.gca().set_xticklabels(plt.gca().get_xticks(),weight='bold')
    plt.gca().set_yticklabels(plt.gca().get_yticks(),weight='bold')
    plt.show()
    txt=str(t+1)
    plt.savefig('/home/awing/orcestra/precipratesnaps/PrecipRate_'+txt.zfill(3)+'.png')
    #plt.savefig('/home/awing/orcestra/OLRsnaps/OLR_'+txt.zfill(3)+'.png')
    plt.close()
