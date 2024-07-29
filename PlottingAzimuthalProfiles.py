#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from datetime import datetime

def az(xord,yord,data,xcen,ycen,r):

	"""
	Function to take, of data about a center location
	azmean(xord,yord,data,xcen,ycen,rr)
	
	Input: 
	xord: x-axis values 
	yord: y-axis values 
	data: data as function (xord,yord) 
	xcen: x value of center
	ycen: y value of center
	r: desired array of radii as 2D array with dimension 1,nr
	
	Output: 
	azmean:, value
	"""
        
	nr = r.shape[1] #size of radius array. right now it's 2D so I am selecting the correct dimension
	nslice = 361 # number of theta slices

    #pre allocate the data
	az = np.zeros((nslice,nr))
	xcircle = np.zeros((nslice,nr))
	ycircle = np.zeros((nslice,nr))
	datai=np.empty((nslice,nr)) #create empty array for the second function
    
    # 360 degrees around the circle 
	dtheta = 360/nslice
	theta = np.atleast_2d(np.linspace(0,360-dtheta,num=nslice)) #array from 0 to 360 with dtheta steps
    # the atleast_2d creates a 2D, with 1,1000 dimension. I am doing that to do the dot product
    # and broadcast the operation.

    #interp2 outputs another function F where you can call xcircle and ycircle
	F = interpolate.interp2d(xord, yord, data) # use the function to interpolate the data 

    # because I have the 2D arrays, and they have one dimension in common (1), I can do 
    # the "dot" product. The important thing is to get the right order and transpose order
    # it took me a couple of tries.
	xcircle = np.array(xcen) + (np.cos(np.deg2rad(theta.transpose()))*r) 
	ycircle = np.array(ycen) + (np.sin(np.deg2rad(theta.transpose()))*r) 

    #When you pass the whole F(xcircle[:,ix],ycircle[:,ix])n
    # you get a 1000x1000 array, and diagonal or firstrow or firstcolumn won't work. so need to loop
	for ix in range(0,nr): 
		for i_nslice in range(0, nslice):
			datai[i_nslice,ix]=F(xcircle[i_nslice,ix],ycircle[i_nslice,ix])
    
    #take azimuthal mean (ignoring cosine of latitude)
	az = datai
    
	return az

def calcshipHeading(shipLon,shipLat):
    dlon = np.diff(shipLon)
    shipHeading=np.empty(len(dlon))
    for hh in range(len(dlon)):
        X = np.cos(np.radians(shipLat[hh+1]))*np.sin(np.radians(dlon[hh]))
        Y = np.cos(np.radians(shipLat[hh]))*np.sin(np.radians(shipLat[hh+1]))-np.sin(np.radians(shipLat[hh]))*np.cos(np.radians(shipLat[hh+1]))*np.cos(np.radians(dlon[hh]))
        init_bearing = np.arctan2(X,Y)
        shipHeading[hh] = np.degrees(init_bearing)

    ideg = np.squeeze(np.where(shipHeading<0))
    shipHeading[ideg]=shipHeading[ideg]+360

    return shipHeading
    
#For purposes of ship heading, calculate the interpolated clons and clats of the ship based on the original track data
#FIELD CAMPAIGN TRACK
trackstarttime = "2016-08-10 12:00:00"
trackendtime = "2016-09-23 20:00:00"
#Set the start time as hours since 1900-01-01 00:00:00
starttime=datetime(1900,1,1,0,0,0,0)
#Get hours since 1900-01-01 00:00:00, then add an hour for each until you get to 2016-09-23 20:00:00
trackstarttime = int(abs(starttime-datetime(2016,8,10,12,0,0,0)).total_seconds() / 3600)
trackendtime = int(abs(starttime-datetime(2016,9,23,20,0,0,0)).total_seconds() / 3600)

hr = 12
d = 10
mo = 8
year = 2016
strTimes = []
for t in range(0,(trackendtime+1 - trackstarttime), 1):
    if(t==0):
        strTimes.append(str(datetime(year,mo,d,hr,0,0,0)))
        hr = hr+1
    else:
        if(hr!=23):
            strTimes.append(str(datetime(year,mo,d,hr,0,0,0)))
            hr = hr+1
        elif(hr==23 and d!=31):
            strTimes.append(str(datetime(year,mo,d,hr,0,0,0)))
            hr = 0
            d = d+1
        elif(hr==23 and d==31):
            strTimes.append(str(datetime(year,mo,d,hr,0,0,0)))
            hr = 0
            d = 1
            mo = mo+1

#Interpolate the lats and lons to length of hourly time array
#Clats
oldclats = [16.88, 14.5, -0.5, 14.5, 2.5, 7.5, 7.5, 14.5, 2.5, 14.5, 7.5, 7.5, 14.5, 5, 7.9, 14.9, 13.15]
oldclats = np.array(oldclats)

x = np.arange(oldclats.size)  # [0, 1, 2, 3, 4, 5]

x_stretch = np.linspace(
    start=x[0], stop=x[-1], num=len(strTimes),)  # [0, 0.625, 1.25, 1.875, 2.5, 3.125, 3.75, 4.375, 5]

oldclats_stretch = np.interp(x_stretch, x, oldclats)
clats = oldclats_stretch
#Clons
oldclons = [-24.98, -23, -23, -23.4, -23.4, -23.8, -32, -32, -32.4, -32.4, -32.8, -47, -47, -47.4, -47.4, -47.4, -59.42]
oldclons = np.array(oldclons)

x = np.arange(oldclons.size)  # [0, 1, 2, 3, 4, 5]

x_stretch = np.linspace(
    start=x[0], stop=x[-1], num=len(strTimes),)  # [0, 0.625, 1.25, 1.875, 2.5, 3.125, 3.75, 4.375, 5]

oldclons_stretch = np.interp(x_stretch, x, oldclons)
clons = oldclons_stretch

#Open data file
data = xr.open_dataset('./ERA-5_CampaignTrack.nc')

#Plot a variable of your choice
#Chose precip and converting from m/s to mm/hr
precipdata = data.TotPrecipratebox * 1000 * 3600
latdata = data.latitude
londata = data.longitude
clatdata = data.centerLat
clondata = data.centerLon
#Calculate the ship heading array
shipheadingdata  = calcshipHeading(clons,clats)

#Loop through all the track points
for t in range(0,len(precipdata)):
    #Theta slices
    thetas = np.radians(np.arange(0, 361, 1))
    #Set radius for azimuth function
    r = np.atleast_2d(np.array([0, 0.25, 0.5, 0.75, 1, 1.25]))  
    nr = r.shape[1]
    radius = np.array([0, 0.25, 0.5, 0.75, 1, 1.25])
    #Set the data variable
    precipbox = precipdata[t]
    #Put the data through the azimuth function
    #Get the clat, clon, latbox, lonbox at the specified time
    clat = clatdata[t]
    clon = clondata[t]
    lats = latdata[t]
    lons = londata[t]
    #To get the correct output from az function, you need to switch the inputs of lat/lon accordingly
    precipaz = az(lats,lons,np.swapaxes(precipbox,0,1),clat,clon,r)
    #Use the ship heading at the time step to determine what to block out
    SH = shipheadingdata[t]
    #Start of blocking piece
    startblock = int(SH + 90)
    endblock = int(SH + 210)
    #Make sure that starting/ending block is within 0-360deg
    if(startblock>360):
        startblock = int(startblock-360)
    if(endblock>360):
        endblock = int(endblock-360)
    #NaN out the blocked area using the start and end blocks
    if(endblock < startblock):
        precipaz[startblock:,:] = precipaz[startblock:,:] * np.nan
        precipaz[0:endblock+1,:] = precipaz[0:endblock+1,:] * np.nan
    else:
        precipaz[startblock:endblock+1,:] = precipaz[startblock:endblock+1,:] * np.nan
    #Start assembling the plotting
    plt.rcParams["font.weight"]="bold"
    plt.rcParams["axes.labelweight"]="bold"
    fig = plt.figure(figsize=[20,10])
    ax = fig.add_subplot(111, polar=True)
    #Assemble the polar plot
    precippolar = ax.contourf(thetas, radius, np.transpose(precipaz), cmap='YlGn', levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,5.5,10.5,15.5,20.5,25.5,30.5,35.5], origin='lower')
    #Set the directions so 0 is N, 90 is E, etc.
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(180)
    #Plot the color bar
    cbar = fig.colorbar(precippolar)
    cbar.set_label("Precip. Rate [mm/hr]", fontweight='bold')
    cbar.ax.set_yticklabels(cbar.ax.get_yticks(),weight='bold')
    #Set the title
    ax.set_title('Precip Rate by Azimuth', fontsize=18, weight='bold', loc='center')
    #Add an arrow to the plot to show the heading direction
    ax.arrow(SH/180.*np.pi,0,0,0.5, facecolor='black',edgecolor='black', width=0.1, linewidth=0.2, alpha=0.5)
    #Save the figure
    txt = str(t+1)
    plt.savefig('./hrlyprecipratesnaps/ERA-5_PrecipRate_AzSnapshot_NoLSmask_'+txt.zfill(3)+'.png')
    #plt.show()
    plt.close()

    #Plot Where the Interpolated track and actual track is
    #Plot the ship tracks
    #Create figure of world map
    plt.figure(figsize=[20,10])
    #Makes everything bold
    plt.rcParams["font.weight"]="bold"
    plt.rcParams["axes.labelweight"]="bold"
    #Make Mercator plot
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    #Plot tracks
    #Interpolated
    ax.plot(clons, clats, linestyle='-', linewidth=3,color='black', transform=ccrs.PlateCarree())
    #Scatter start and end points
    ax.scatter(clons[0],clats[0], marker='x', s=500,color='red', transform=ccrs.PlateCarree())
    ax.scatter(clons[-1],clats[-1], marker='o', s=500,color='red', facecolor='white', transform=ccrs.PlateCarree())
    #Scatter a given point along the track
    ax.scatter(clons[t],clats[t], marker='o', s=250,color='red', transform=ccrs.PlateCarree())
    #Set the lat/lon axes
    ax.set_xticks([270,285,300,315,330,345,359.999999999], crs=ccrs.PlateCarree())
    ax.set_yticks([-30,-15,0,15,30], crs=ccrs.PlateCarree())
    #Make tick labels larger
    ax.tick_params(axis='both',labelsize=15)
    #Have plot only go from 30S to 30N
    ax.set_ylim(-30,30)
    #Format the Lat/Lon correctly for respective axes
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.xaxis.set_major_formatter(lon_formatter)
    #Put yaxis label on RHS
    ax.yaxis.set_ticks_position('right')
    #Grid the lats/lons
    ax.grid(linestyle='--',color='lightgrey',linewidth=1.5, alpha=0.75)
    #Add Coastlines
    ax.coastlines(resolution='110m',linewidth=2)
    plt.title('Field Campaign Interpolated Ship Track\nX is the Start, O is the Finish', fontweight='bold', fontsize=20)
    #Save and show plot
    txt = str(t+1)
    plt.savefig('./ShipTrack/FieldCampaign_ShipLocationAlongtrack_'+txt.zfill(3)+'.png')
    #plt.show()
    plt.close()