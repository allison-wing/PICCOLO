#%%
from re import S
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp2d
import scipy.stats as stats
from datetime import datetime

########################MATH FUNCTIONS(box averaging, etc.)############################################
def boxavg(thing,lat,lon):
    coslat_values = np.transpose(np.tile(np.cos(np.deg2rad(lat)),(len(lon),1)))
    thing1 = thing*coslat_values
    thing2 = thing1/thing1
    average = np.nansum(np.nansum(thing1,0))/np.nansum(np.nansum(coslat_values*thing2,0))

    return average

def azmean(xord,yord,data,xcen,ycen,r):

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
	nslice = 1000 # number of theta slices

    #pre allocate the data
	azmean = np.zeros(nr)
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
	azmean = np.nanmean(datai, axis=0)
    
	return azmean

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
	nslice = 1000 # number of theta slices

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

#File directories
era5MainDir = '/mars/tank3/era5/'
RadFluxDir = 'forecasts/RAD/redownloaded_fluxes/' #files are by month ex) 2007-06.nc
SfcFluxDir = 'forecasts/RAD/lhf_shf/' #files are by month ex) 2007-06.nc
PrecipDir = 'forecasts/precip/' #files are by month ex) 2007-06.nc, hrs since 01-01-1900
SfcVarsDir = 'surface_variables/' #files are by month ex) 2007-06.nc
WindboxDir = 'uv_10m_wind/' #files are by month ex) 2007-06.nc
vortdir = '/mars/tank3/era5/pressure_variables/relative_vorticity/' #One file every 6hrs (ex. rv_850.nc, var is 'vo' s^-1, atmosphere relative vort) hrs since 01-01-1900
geopotdir = '/mars/tank3/era5/pressure_variables/geopotential/' #One file every 6hrs, 500hpa only (geopotential_500.nc, var is 'z' m^2/s^2) hrs since 01-01-1900

#Open rel. vort files here as they have all the times in them
Vort900_file = xr.open_dataset(vortdir+'rv_900.nc', decode_times=False, use_cftime=False)
Vort850_file = xr.open_dataset(vortdir+'rv_850.nc', decode_times=False, use_cftime=False)
Vort600_file = xr.open_dataset(vortdir+'rv_600.nc', decode_times=False, use_cftime=False)
#All have same times so just pick one to get t indexing lists
vorttarr = Vort850_file.time
vortitarr = pd.Index(vorttarr)
vortitlst = vortitarr.tolist()

#Open the lsm and grab the variable
lsmask = xr.open_dataset('/mars/tank3/era5/landmask/land-sea_mask.nc')
#Now gather and put general lats/lons list into index format to use later for gathering lat/lonbox data for a given time
lats = np.array(lsmask['latitude'])
lons = np.array(lsmask['longitude'])
ilats = pd.Index(lats)
ilons = pd.Index(lons)
ilatlist = ilats.tolist()
ilonlist = ilons.tolist()
lsm = lsmask.lsm.isel(time=0)
lsmask.close()

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
clons =clons + 360

#Get variable lengths set up
latres = 0.25
lonres = 0.25
r = np.atleast_2d(np.array([0, 0.25, 0.5, 0.75, 1, 1.25]))  
nr = r.shape[1]
nslice = 1000

#Latitude, Longitude amounts to get 1.5X1.5 deg box
latlen = int(2.5/latres + 1) #Center lat position is one index, then 0.75 degrees up and 0.75 degrees down
lonlen = int(2.5/lonres + 1) #Center lon position is one index, then 0.75 degrees left and 0.75 degrees right
numsteps = len(strTimes)
#Get all the variable save names
#Create the 4-D arrays for all the variables desired
windbox_azsave = np.ones((numsteps, nslice, nr))*np.nan

sfcpres_azsave = np.ones((numsteps, nslice, nr))*np.nan

totprecip_azsave = np.ones((numsteps, nslice, nr))*np.nan

RelVort900_azsave = np.ones((numsteps, nslice, nr))*np.nan

RelVort850_azsave = np.ones((numsteps, nslice, nr))*np.nan

RelVort600_azsave = np.ones((numsteps, nslice, nr))*np.nan

#GeoPot500_save = np.ones((numstorms, numsteps, latlen, lonlen))*np.nan

NetLWSfc_azsave = np.ones((numsteps, nslice, nr))*np.nan

NetLWSfcCS_azsave = np.ones((numsteps, nslice, nr))*np.nan

NetLWToa_azsave = np.ones((numsteps, nslice, nr))*np.nan

NetLWToaCS_azsave = np.ones((numsteps, nslice, nr))*np.nan

ClmnLWfluxConv_azsave = np.ones((numsteps, nslice, nr))*np.nan

NetSWSfc_azsave = np.ones((numsteps, nslice, nr))*np.nan

NetSWSfcCS_azsave = np.ones((numsteps, nslice, nr))*np.nan

NetSWToa_azsave = np.ones((numsteps, nslice, nr))*np.nan

NetSWToaCS_azsave = np.ones((numsteps, nslice, nr))*np.nan

ClmnSWfluxConv_azsave = np.ones((numsteps, nslice, nr))*np.nan

ClmnLWfluxConvCS_azsave = np.ones((numsteps, nslice, nr))*np.nan

ClmnSWfluxConvCS_azsave = np.ones((numsteps, nslice, nr))*np.nan

ClmnRadfluxConv_azsave = np.ones((numsteps, nslice, nr))*np.nan

ClmnRadfluxConvCS_azsave = np.ones((numsteps, nslice, nr))*np.nan

OLR_azsave = np.ones((numsteps, nslice, nr))*np.nan

hfls_azsave = np.ones((numsteps, nslice, nr))*np.nan

hfss_azsave = np.ones((numsteps, nslice, nr))*np.nan

sfcMoistEnthalpyFlux_azsave = np.ones((numsteps, nslice, nr))*np.nan

#Non-Radiative variables (ex: slp, wind, years, months, days, hours, latbox, lonbox, clat, clon, etc.)
#3D variables
h_azmean = np.ones((numsteps,nr))*np.nan

LHF_azmean = np.ones((numsteps,nr))*np.nan

SHF_azmean = np.ones((numsteps,nr))*np.nan

#CpTdiff_azmean = np.ones((nstorms,ntracks,nr))*np.nan

#Tdiff_azmean = np.ones((nstorms,ntracks,nr))*np.nan

#LvQdiff_azmean = np.ones((nstorms,ntracks,nr))*np.nan

#Qdiff_azmean = np.ones((nstorms,ntracks,nr))*np.nan

wind_azmean = np.ones((numsteps,nr))*np.nan

totprecip_azmean = np.ones((numsteps,nr))*np.nan

windbin_azmean = np.ones((numsteps,nr))*np.nan

#2D variables
LW_SEF_ratio = np.ones((numsteps))*np.nan

maxwind_save = np.ones((numsteps))*np.nan

minSLP_save = np.ones((numsteps))*np.nan

Clat_save = np.ones((numsteps))*np.nan

Clon_save = np.ones((numsteps))*np.nan

year_save = np.ones((numsteps))*np.nan

month_save = np.ones((numsteps))*np.nan

day_save = np.ones((numsteps))*np.nan

hour_save = np.ones((numsteps))*np.nan
#%%
for t, time in enumerate(strTimes):
    print(t)
    Mostr = time[5:7]
    mo = int(Mostr)
    Dstr = time[8:10]
    d = int(Dstr)
    Hrstr = time[11:13]
    hr = int(Hrstr)
    #Open the surface var and h data file
    SfcVarFile = xr.open_dataset(era5MainDir+SfcVarsDir+str(year)+'-'+Mostr+'.nc', decode_times=False, use_cftime=False)

    #Get the specific vars and then close the parent dataset after
    #Get the clat/clon position that is closest to what is provided in track data
    clat = SfcVarFile['latitude'].sel(latitude=clats[t], method='nearest')
    clon = SfcVarFile['longitude'].sel(longitude=clons[t], method='nearest')
    #Get the index of the above found clat/clon
    iclat = ilatlist.index(clat)
    iclon = ilonlist.index(clon)
    #Now set up bounds of 10X10 deg box based on index spacing and must go 1 higher for largest bound (For GFDL-CM4 it is 1deg lat X 1.25deg lon so +/-5 and +/-4 respectively)
    latmax = iclat+int((latlen-1)/2+1)
    latmin = iclat-int((latlen-1)/2)
    lonmax = iclon+int((lonlen-1)/2+1)
    lonmin = iclon-int((lonlen-1)/2)
    #Now gather the lat/lon array for the box
    latbox = np.array(SfcVarFile.latitude.isel(latitude=slice(latmin,latmax)))
    #Now set up if there is a lon bound issue (right or left)
    #Right lon bound issue
    if(iclon>1419):
        excessright = lonmax-1440
        rightlonmin = 0
        rightlonmax = excessright
        leftlonmin = lonmin
        leftlonmax = lonmax
        rightlonbox = np.array(SfcVarFile.longitude.isel(longitude=slice(rightlonmin,rightlonmax)))
        leftlonbox = np.array(SfcVarFile.longitude.isel(longitude=slice(lonmin,lonmax)))
        lonbox = np.concatenate([leftlonbox,rightlonbox],axis=0)
    #Left lon bound issue
    elif(iclon<20):
        excessleft = np.abs(lonmin)
        rightlonmin=0
        rightlonmax = lonmax
        leftlonmin = 1440-excessleft
        leftlonmax = 1440
        leftlonbox = np.array(SfcVarFile.longitude.isel(longitude=slice(leftlonmin,leftlonmax)))
        rightlonbox = np.array(SfcVarFile.longitude.isel(longitude=slice(rightlonmin,rightlonmax)))
        lonbox = np.concatenate([leftlonbox,rightlonbox],axis=0)
    #No lon bound issue
    else:
        #Will just make the right lonmin/max = 0 and have it be an empty array
        leftlonmin=lonmin
        leftlonmax=lonmax
        rightlonmin=0
        rightlonmax=0
        lonbox = np.array(SfcVarFile.longitude.isel(longitude=slice(leftlonmin,leftlonmax)))
    #Now make a 2D array based on the land-sea mask that is zeros or NaN if >20%
    landsea_zerosNaNs = np.zeros((len(latbox),len(lonbox)))
    #Open the parent land-sea mask file from outside the loop that is sliced according to the lat/lon bounds above
    landsea_slicedleft = np.array(lsm.isel(latitude=slice(latmin,latmax),longitude=slice(leftlonmin,leftlonmax)))
    landsea_slicedright = np.array(lsm.isel(latitude=slice(latmin,latmax),longitude=slice(rightlonmin,rightlonmax)))
    landsea_sliced = np.concatenate([landsea_slicedleft,landsea_slicedright],axis=1)
    #Now loop through the sliced land-sea mask to assign NaNs to the grid points that are >20
    for i in range(0,len(latbox)):
        for j in range(0,len(lonbox)):
            if(landsea_sliced[i][j] > 0.2):
                landsea_zerosNaNs[i][j]=np.nan
    
    landsea_zerosNaNsaz = az(lonbox,latbox,landsea_zerosNaNs,clon,clat,r)
    
    #Now get the endtime1 and endtime2 by case depending on hour 
    #00z case for all other months(need the last index of the prior month's file (23z) and the second index of the current month's file (01z))
    if(hr==0 and d==1):
        #Get the prior month for the first file
        mo1 = mo-1
        if(mo1<10):
            mo1str = '0'+str(mo1)
        else:
            mo1str = str(mo1)
        #Get the current month for the first file
        mo2 = mo
        if(mo2<10):
            mo2str = '0'+str(mo2)
        else:
            mo2str = str(mo2)
        filetime1 = str(year)+'-'+mo1str+'.nc'
        filetime2 = str(year)+'-'+mo2str+'.nc'

        #Open the flux files
        SfcFluxFile1 = xr.open_dataset(era5MainDir+SfcFluxDir+filetime1, decode_times=False, use_cftime=False)
        RadFluxFile1 = xr.open_dataset(era5MainDir+RadFluxDir+filetime1, decode_times=False, use_cftime=False)
        SfcFluxFile2 = xr.open_dataset(era5MainDir+SfcFluxDir+filetime2, decode_times=False, use_cftime=False)
        RadFluxFile2 = xr.open_dataset(era5MainDir+RadFluxDir+filetime2, decode_times=False, use_cftime=False)
        #Get the times used for the time indexing
        t1 = int(abs(starttime-datetime(year,mo1,31,23,0,0,0)).total_seconds() / 3600)
        t2 = int(abs(starttime-datetime(year,mo2,1,1,0,0,0)).total_seconds() / 3600)
        #Times for file 1
        tarray1 = SfcFluxFile1.time
        itarray1 = pd.Index(tarray1)
        itlist1 = itarray1.tolist()
        #Times for file 2
        tarray2 = SfcFluxFile2.time
        itarray2 = pd.Index(tarray2)
        itlist2 = itarray2.tolist()
        #The indexing
        tind1 = itlist1.index(t1)
        tind2 = itlist2.index(t2)
        #Fluxes
        #hfls
        hfls1left = np.array(SfcFluxFile1.slhf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfls1right = np.array(SfcFluxFile1.slhf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfls1 =  np.concatenate([hfls1left,hfls1right],axis=1) 
        hfls2left = np.array(SfcFluxFile2.slhf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfls2right = np.array(SfcFluxFile2.slhf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfls2 =  np.concatenate([hfls2left,hfls2right],axis=1)
        hfls = (hfls1 + hfls2) / (-3600*2)
        #hfss
        hfss1left = np.array(SfcFluxFile1.sshf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfss1right = np.array(SfcFluxFile1.sshf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfss1 =  np.concatenate([hfss1left,hfss1right],axis=1) 
        hfss2left = np.array(SfcFluxFile2.sshf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfss2right = np.array(SfcFluxFile2.sshf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfss2 =  np.concatenate([hfss2left,hfss2right],axis=1)
        hfss = (hfss1 + hfss2) / (-3600*2)
        #SW Sfc
        NetSWsfcleft1 = np.array(RadFluxFile1.ssr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcright1 = np.array(RadFluxFile1.ssr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfc1 = np.concatenate([NetSWsfcleft1,NetSWsfcright1],axis=1)
        NetSWsfcleft2 = np.array(RadFluxFile2.ssr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcright2 = np.array(RadFluxFile2.ssr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfc2 = np.concatenate([NetSWsfcleft2,NetSWsfcright2],axis=1)
        NetSWsfc = (NetSWsfc1 + NetSWsfc2)/(3600*2)
        #SW CS Sfc
        NetSWsfcCSleft1 = np.array(RadFluxFile1.ssrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcCSright1 = np.array(RadFluxFile1.ssrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfcCS1 = np.concatenate([NetSWsfcCSleft1,NetSWsfcCSright1],axis=1)
        NetSWsfcCSleft2 = np.array(RadFluxFile2.ssrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcCSright2 = np.array(RadFluxFile2.ssrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfcCS2 = np.concatenate([NetSWsfcCSleft2,NetSWsfcCSright2],axis=1)
        NetSWsfcCS = (NetSWsfcCS1 + NetSWsfcCS2)/(3600*2)
        #LW Sfc
        NetLWsfcleft1 = np.array(RadFluxFile1.str.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcright1 = np.array(RadFluxFile1.str.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfc1 = np.concatenate([NetLWsfcleft1,NetLWsfcright1],axis=1)
        NetLWsfcleft2 = np.array(RadFluxFile2.str.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcright2 = np.array(RadFluxFile2.str.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfc2 = np.concatenate([NetLWsfcleft2,NetLWsfcright2],axis=1)
        NetLWsfc = (NetLWsfc1 + NetLWsfc2)/(3600*2)
        #LW CS Sfc
        NetLWsfcCSleft1 = np.array(RadFluxFile1.strc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcCSright1 = np.array(RadFluxFile1.strc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfcCS1 = np.concatenate([NetLWsfcCSleft1,NetLWsfcCSright1],axis=1)
        NetLWsfcCSleft2 = np.array(RadFluxFile2.strc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcCSright2 = np.array(RadFluxFile2.strc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfcCS2 = np.concatenate([NetLWsfcCSleft2,NetLWsfcCSright2],axis=1)
        NetLWsfcCS = (NetLWsfcCS1 + NetLWsfcCS2)/(3600*2)
        #SW TOA
        NetSWToaleft1 = np.array(RadFluxFile1.tsr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaright1 = np.array(RadFluxFile1.tsr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToa1 = np.concatenate([NetSWToaleft1,NetSWToaright1],axis=1)
        NetSWToaleft2 = np.array(RadFluxFile2.tsr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaright2 = np.array(RadFluxFile2.tsr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToa2 = np.concatenate([NetSWToaleft2,NetSWToaright2],axis=1)
        NetSWToa = (NetSWToa1 + NetSWToa2)/(3600*2)
        #SW CS TOA
        NetSWToaCSleft1 = np.array(RadFluxFile1.tsrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaCSright1 = np.array(RadFluxFile1.tsrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToaCS1 = np.concatenate([NetSWToaCSleft1,NetSWToaCSright1],axis=1)
        NetSWToaCSleft2 = np.array(RadFluxFile2.tsrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaCSright2 = np.array(RadFluxFile2.tsrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToaCS2 = np.concatenate([NetSWToaCSleft2,NetSWToaCSright2],axis=1)
        NetSWToaCS = (NetSWToaCS1 + NetSWToaCS2)/(3600*2)
        #LW TOA
        NetLWToaleft1 = np.array(RadFluxFile1.ttr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaright1 = np.array(RadFluxFile1.ttr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToa1 = np.concatenate([NetLWToaleft1,NetLWToaright1],axis=1)
        NetLWToaleft2 = np.array(RadFluxFile2.ttr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaright2 = np.array(RadFluxFile2.ttr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToa2 = np.concatenate([NetLWToaleft2,NetLWToaright2],axis=1)
        NetLWToa = (NetLWToa1 + NetLWToa2)/(3600*2)
        #LW CS TOA
        NetLWToaCSleft1 = np.array(RadFluxFile1.ttrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaCSright1 = np.array(RadFluxFile1.ttrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToaCS1 = np.concatenate([NetLWToaCSleft1,NetLWToaCSright1],axis=1)
        NetLWToaCSleft2 = np.array(RadFluxFile2.ttrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaCSright2 = np.array(RadFluxFile2.ttrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToaCS2 = np.concatenate([NetLWToaCSleft2,NetLWToaCSright2],axis=1)
        NetLWToaCS = (NetLWToaCS1 + NetLWToaCS2)/(3600*2)
        #Close the flux files
        RadFluxFile1.close()
        SfcFluxFile1.close()
        RadFluxFile2.close()
        SfcFluxFile2.close()

        hflsaz = az(lonbox,latbox,hfls,clon,clat,r)
        hfssaz = az(lonbox,latbox,hfss,clon,clat,r)
        NetSWsfcaz = az(lonbox,latbox,NetSWsfc,clon,clat,r)
        NetSWsfcCSaz = az(lonbox,latbox,NetSWsfcCS,clon,clat,r)
        NetLWsfcaz = az(lonbox,latbox,NetLWsfc,clon,clat,r)
        NetLWsfcCSaz = az(lonbox,latbox,NetLWsfcCS,clon,clat,r)
        NetSWToaaz = az(lonbox,latbox,NetSWToa,clon,clat,r)
        NetSWToaCSaz = az(lonbox,latbox,NetSWToaCS,clon,clat,r)
        NetLWToaaz = az(lonbox,latbox,NetLWToa,clon,clat,r)
        NetLWToaCSaz = az(lonbox,latbox,NetLWToaCS,clon,clat,r)
    #23z case for end of month
    elif(hr==23 and d==31):
        #Get the prior month for the first file
        mo1 = mo
        if(mo1<10):
            mo1str = '0'+str(mo1)
        else:
            mo1str = str(mo1)
        #Get the current month for the first file
        mo2 = mo + 1
        if(mo2<10):
            mo2str = '0'+str(mo2)
        else:
            mo2str = str(mo2)
        filetime1 = str(year)+'-'+mo1str+'.nc'
        filetime2 = str(year)+'-'+mo2str+'.nc'

        #Open the flux files
        SfcFluxFile1 = xr.open_dataset(era5MainDir+SfcFluxDir+filetime1, decode_times=False, use_cftime=False)
        RadFluxFile1 = xr.open_dataset(era5MainDir+RadFluxDir+filetime1, decode_times=False, use_cftime=False)
        SfcFluxFile2 = xr.open_dataset(era5MainDir+SfcFluxDir+filetime2, decode_times=False, use_cftime=False)
        RadFluxFile2 = xr.open_dataset(era5MainDir+RadFluxDir+filetime2, decode_times=False, use_cftime=False)
        #Get the times used for the time indexing
        t1 = int(abs(starttime-datetime(year,mo1,31,22,0,0,0)).total_seconds() / 3600)
        t2 = int(abs(starttime-datetime(year,mo2,1,0,0,0,0)).total_seconds() / 3600)
        #Times for file 1
        tarray1 = SfcFluxFile1.time
        itarray1 = pd.Index(tarray1)
        itlist1 = itarray1.tolist()
        #Times for file 2
        tarray2 = SfcFluxFile2.time
        itarray2 = pd.Index(tarray2)
        itlist2 = itarray2.tolist()
        #The indexing
        tind1 = itlist1.index(t1)
        tind2 = itlist2.index(t2)
        #Fluxes
        #hfls
        hfls1left = np.array(SfcFluxFile1.slhf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfls1right = np.array(SfcFluxFile1.slhf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfls1 =  np.concatenate([hfls1left,hfls1right],axis=1) 
        hfls2left = np.array(SfcFluxFile2.slhf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfls2right = np.array(SfcFluxFile2.slhf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfls2 =  np.concatenate([hfls2left,hfls2right],axis=1)
        hfls = (hfls1 + hfls2) / (-3600*2)
        #hfss
        hfss1left = np.array(SfcFluxFile1.sshf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfss1right = np.array(SfcFluxFile1.sshf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfss1 =  np.concatenate([hfss1left,hfss1right],axis=1) 
        hfss2left = np.array(SfcFluxFile2.sshf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfss2right = np.array(SfcFluxFile2.sshf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfss2 =  np.concatenate([hfss2left,hfss2right],axis=1)
        hfss = (hfss1 + hfss2) / (-3600*2)
        #SW Sfc
        NetSWsfcleft1 = np.array(RadFluxFile1.ssr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcright1 = np.array(RadFluxFile1.ssr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfc1 = np.concatenate([NetSWsfcleft1,NetSWsfcright1],axis=1)
        NetSWsfcleft2 = np.array(RadFluxFile2.ssr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcright2 = np.array(RadFluxFile2.ssr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfc2 = np.concatenate([NetSWsfcleft2,NetSWsfcright2],axis=1)
        NetSWsfc = (NetSWsfc1 + NetSWsfc2)/(3600*2)
        #SW CS Sfc
        NetSWsfcCSleft1 = np.array(RadFluxFile1.ssrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcCSright1 = np.array(RadFluxFile1.ssrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfcCS1 = np.concatenate([NetSWsfcCSleft1,NetSWsfcCSright1],axis=1)
        NetSWsfcCSleft2 = np.array(RadFluxFile2.ssrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcCSright2 = np.array(RadFluxFile2.ssrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfcCS2 = np.concatenate([NetSWsfcCSleft2,NetSWsfcCSright2],axis=1)
        NetSWsfcCS = (NetSWsfcCS1 + NetSWsfcCS2)/(3600*2)
        #LW Sfc
        NetLWsfcleft1 = np.array(RadFluxFile1.str.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcright1 = np.array(RadFluxFile1.str.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfc1 = np.concatenate([NetLWsfcleft1,NetLWsfcright1],axis=1)
        NetLWsfcleft2 = np.array(RadFluxFile2.str.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcright2 = np.array(RadFluxFile2.str.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfc2 = np.concatenate([NetLWsfcleft2,NetLWsfcright2],axis=1)
        NetLWsfc = (NetLWsfc1 + NetLWsfc2)/(3600*2)
        #LW CS Sfc
        NetLWsfcCSleft1 = np.array(RadFluxFile1.strc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcCSright1 = np.array(RadFluxFile1.strc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfcCS1 = np.concatenate([NetLWsfcCSleft1,NetLWsfcCSright1],axis=1)
        NetLWsfcCSleft2 = np.array(RadFluxFile2.strc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcCSright2 = np.array(RadFluxFile2.strc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfcCS2 = np.concatenate([NetLWsfcCSleft2,NetLWsfcCSright2],axis=1)
        NetLWsfcCS = (NetLWsfcCS1 + NetLWsfcCS2)/(3600*2)
        #SW TOA
        NetSWToaleft1 = np.array(RadFluxFile1.tsr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaright1 = np.array(RadFluxFile1.tsr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToa1 = np.concatenate([NetSWToaleft1,NetSWToaright1],axis=1)
        NetSWToaleft2 = np.array(RadFluxFile2.tsr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaright2 = np.array(RadFluxFile2.tsr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToa2 = np.concatenate([NetSWToaleft2,NetSWToaright2],axis=1)
        NetSWToa = (NetSWToa1 + NetSWToa2)/(3600*2)
        #SW CS TOA
        NetSWToaCSleft1 = np.array(RadFluxFile1.tsrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaCSright1 = np.array(RadFluxFile1.tsrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToaCS1 = np.concatenate([NetSWToaCSleft1,NetSWToaCSright1],axis=1)
        NetSWToaCSleft2 = np.array(RadFluxFile2.tsrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaCSright2 = np.array(RadFluxFile2.tsrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToaCS2 = np.concatenate([NetSWToaCSleft2,NetSWToaCSright2],axis=1)
        NetSWToaCS = (NetSWToaCS1 + NetSWToaCS2)/(3600*2)
        #LW TOA
        NetLWToaleft1 = np.array(RadFluxFile1.ttr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaright1 = np.array(RadFluxFile1.ttr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToa1 = np.concatenate([NetLWToaleft1,NetLWToaright1],axis=1)
        NetLWToaleft2 = np.array(RadFluxFile2.ttr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaright2 = np.array(RadFluxFile2.ttr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToa2 = np.concatenate([NetLWToaleft2,NetLWToaright2],axis=1)
        NetLWToa = (NetLWToa1 + NetLWToa2)/(3600*2)
        #LW CS TOA
        NetLWToaCSleft1 = np.array(RadFluxFile1.ttrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaCSright1 = np.array(RadFluxFile1.ttrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToaCS1 = np.concatenate([NetLWToaCSleft1,NetLWToaCSright1],axis=1)
        NetLWToaCSleft2 = np.array(RadFluxFile2.ttrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaCSright2 = np.array(RadFluxFile2.ttrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToaCS2 = np.concatenate([NetLWToaCSleft2,NetLWToaCSright2],axis=1)
        NetLWToaCS = (NetLWToaCS1 + NetLWToaCS2)/(3600*2)
        #Close the flux files
        RadFluxFile1.close()
        SfcFluxFile1.close()
        RadFluxFile2.close()
        SfcFluxFile2.close()

        hflsaz = az(lonbox,latbox,hfls,clon,clat,r)
        hfssaz = az(lonbox,latbox,hfss,clon,clat,r)
        NetSWsfcaz = az(lonbox,latbox,NetSWsfc,clon,clat,r)
        NetSWsfcCSaz = az(lonbox,latbox,NetSWsfcCS,clon,clat,r)
        NetLWsfcaz = az(lonbox,latbox,NetLWsfc,clon,clat,r)
        NetLWsfcCSaz = az(lonbox,latbox,NetLWsfcCS,clon,clat,r)
        NetSWToaaz = az(lonbox,latbox,NetSWToa,clon,clat,r)
        NetSWToaCSaz = az(lonbox,latbox,NetSWToaCS,clon,clat,r)
        NetLWToaaz = az(lonbox,latbox,NetLWToa,clon,clat,r)
        NetLWToaCSaz = az(lonbox,latbox,NetLWToaCS,clon,clat,r)

    #00z (not on day 1 of month), and all other hours
    else:
        #Check the hour situation
        if(hr==0):
            hr1=23
            hr2 = hr+1
            d1=d-1
            d2 = d
        elif(hr==23):
            hr1=hr-1
            hr2 = 0
            d1=d
            d2 = d + 1
        elif(hr!=23 and hr!=0):
            hr1 = hr - 1
            hr2 = hr + 1
            d1 = d
            d2 = d
        if(mo<10):
            mostr='0'+str(mo)
        else:
            mostr=str(mo)
        filetime = str(year)+'-'+mostr+'.nc'
        SfcFluxFile = xr.open_dataset(era5MainDir+SfcFluxDir+filetime, decode_times=False, use_cftime=False)
        RadFluxFile = xr.open_dataset(era5MainDir+RadFluxDir+filetime, decode_times=False, use_cftime=False)
        #Get the times used for the time indexing
        t1 = int(abs(starttime-datetime(year,mo,d1,hr1,0,0,0)).total_seconds() / 3600)
        t2 = int(abs(starttime-datetime(year,mo,d2,hr2,0,0,0)).total_seconds() / 3600)
        #Times for file
        tarray = SfcFluxFile.time
        itarray = pd.Index(tarray)
        itlist = itarray.tolist()
        #The indexing
        tind1 = itlist.index(t1)
        tind2 = itlist.index(t2)
        #Fluxes
        #hfls
        hfls1left = np.array(SfcFluxFile.slhf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfls1right = np.array(SfcFluxFile.slhf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfls1 =  np.concatenate([hfls1left,hfls1right],axis=1) 
        hfls2left = np.array(SfcFluxFile.slhf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfls2right = np.array(SfcFluxFile.slhf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfls2 =  np.concatenate([hfls2left,hfls2right],axis=1)
        hfls = (hfls1 + hfls2) / (-3600*2)
        hflsaz = az(lonbox,latbox,hfls,clon,clat,r)
        #hfss
        hfss1left = np.array(SfcFluxFile.sshf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfss1right = np.array(SfcFluxFile.sshf.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfss1 =  np.concatenate([hfss1left,hfss1right],axis=1) 
        hfss2left = np.array(SfcFluxFile.sshf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        hfss2right = np.array(SfcFluxFile.sshf.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        hfss2 =  np.concatenate([hfss2left,hfss2right],axis=1)
        hfss = (hfss1 + hfss2) / (-3600*2)
        hfssaz = az(lonbox,latbox,hfss,clon,clat,r)
        #SW Sfc
        NetSWsfcleft1 = np.array(RadFluxFile.ssr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcright1 = np.array(RadFluxFile.ssr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfc1 = np.concatenate([NetSWsfcleft1,NetSWsfcright1],axis=1)
        NetSWsfcleft2 = np.array(RadFluxFile.ssr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcright2 = np.array(RadFluxFile.ssr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfc2 = np.concatenate([NetSWsfcleft2,NetSWsfcright2],axis=1)
        NetSWsfc = (NetSWsfc1 + NetSWsfc2)/(3600*2)
        NetSWsfcaz = az(lonbox,latbox,NetSWsfc,clon,clat,r)
        #SW CS Sfc
        NetSWsfcCSleft1 = np.array(RadFluxFile.ssrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcCSright1 = np.array(RadFluxFile.ssrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfcCS1 = np.concatenate([NetSWsfcCSleft1,NetSWsfcCSright1],axis=1)
        NetSWsfcCSleft2 = np.array(RadFluxFile.ssrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWsfcCSright2 = np.array(RadFluxFile.ssrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWsfcCS2 = np.concatenate([NetSWsfcCSleft2,NetSWsfcCSright2],axis=1)
        NetSWsfcCS = (NetSWsfcCS1 + NetSWsfcCS2)/(3600*2)
        NetSWsfcCSaz = az(lonbox,latbox,NetSWsfcCS,clon,clat,r)
        #LW Sfc
        NetLWsfcleft1 = np.array(RadFluxFile.str.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcright1 = np.array(RadFluxFile.str.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfc1 = np.concatenate([NetLWsfcleft1,NetLWsfcright1],axis=1)
        NetLWsfcleft2 = np.array(RadFluxFile.str.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcright2 = np.array(RadFluxFile.str.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfc2 = np.concatenate([NetLWsfcleft2,NetLWsfcright2],axis=1)
        NetLWsfc = (NetLWsfc1 + NetLWsfc2)/(3600*2)
        NetLWsfcaz = az(lonbox,latbox,NetLWsfc,clon,clat,r)
        #LW CS Sfc
        NetLWsfcCSleft1 = np.array(RadFluxFile.strc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcCSright1 = np.array(RadFluxFile.strc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfcCS1 = np.concatenate([NetLWsfcCSleft1,NetLWsfcCSright1],axis=1)
        NetLWsfcCSleft2 = np.array(RadFluxFile.strc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWsfcCSright2 = np.array(RadFluxFile.strc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWsfcCS2 = np.concatenate([NetLWsfcCSleft2,NetLWsfcCSright2],axis=1)
        NetLWsfcCS = (NetLWsfcCS1 + NetLWsfcCS2)/(3600*2)
        NetLWsfcCSaz = az(lonbox,latbox,NetLWsfcCS,clon,clat,r)
        #SW TOA
        NetSWToaleft1 = np.array(RadFluxFile.tsr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaright1 = np.array(RadFluxFile.tsr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToa1 = np.concatenate([NetSWToaleft1,NetSWToaright1],axis=1)
        NetSWToaleft2 = np.array(RadFluxFile.tsr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaright2 = np.array(RadFluxFile.tsr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToa2 = np.concatenate([NetSWToaleft2,NetSWToaright2],axis=1)
        NetSWToa = (NetSWToa1 + NetSWToa2)/(3600*2)
        NetSWToaaz = az(lonbox,latbox,NetSWToa,clon,clat,r)
        #SW CS TOA
        NetSWToaCSleft1 = np.array(RadFluxFile.tsrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaCSright1 = np.array(RadFluxFile.tsrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToaCS1 = np.concatenate([NetSWToaCSleft1,NetSWToaCSright1],axis=1)
        NetSWToaCSleft2 = np.array(RadFluxFile.tsrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetSWToaCSright2 = np.array(RadFluxFile.tsrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetSWToaCS2 = np.concatenate([NetSWToaCSleft2,NetSWToaCSright2],axis=1)
        NetSWToaCS = (NetSWToaCS1 + NetSWToaCS2)/(3600*2)
        NetSWToaCSaz = az(lonbox,latbox,NetSWToaCS,clon,clat,r)
        #LW TOA
        NetLWToaleft1 = np.array(RadFluxFile.ttr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaright1 = np.array(RadFluxFile.ttr.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToa1 = np.concatenate([NetLWToaleft1,NetLWToaright1],axis=1)
        NetLWToaleft2 = np.array(RadFluxFile.ttr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaright2 = np.array(RadFluxFile.ttr.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToa2 = np.concatenate([NetLWToaleft2,NetLWToaright2],axis=1)
        NetLWToa = (NetLWToa1 + NetLWToa2)/(3600*2)
        NetLWToaaz = az(lonbox,latbox,NetLWToa,clon,clat,r)
        #LW CS TOA
        NetLWToaCSleft1 = np.array(RadFluxFile.ttrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaCSright1 = np.array(RadFluxFile.ttrc.isel(time = tind1, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToaCS1 = np.concatenate([NetLWToaCSleft1,NetLWToaCSright1],axis=1)
        NetLWToaCSleft2 = np.array(RadFluxFile.ttrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(leftlonmin,leftlonmax)))
        NetLWToaCSright2 = np.array(RadFluxFile.ttrc.isel(time = tind2, latitude = slice(latmin,latmax), longitude = slice(rightlonmin,rightlonmax)))
        NetLWToaCS2 = np.concatenate([NetLWToaCSleft2,NetLWToaCSright2],axis=1)
        NetLWToaCS = (NetLWToaCS1 + NetLWToaCS2)/(3600*2)
        NetLWToaCSaz = az(lonbox,latbox,NetLWToaCS,clon,clat,r)
        #Close the flux files
        RadFluxFile.close()
        SfcFluxFile.close()

    #Get the budget set up
    ClmnLWfluxConv = NetLWToa - NetLWsfc
    ClmnLWfluxConvCS = NetLWToaCS - NetLWsfcCS
    ClmnSWfluxConv = NetSWToa - NetSWsfc
    ClmnSWfluxConvCS = NetSWToaCS - NetSWsfcCS
    ClmnRadfluxConv = ClmnLWfluxConv + ClmnSWfluxConv
    ClmnRadfluxConvCS = ClmnLWfluxConvCS + ClmnSWfluxConvCS
    SfcMoistEnthalpyFlux = hfls + hfss
    ClmnLWfluxConvaz = az(lonbox,latbox,ClmnLWfluxConv,clon,clat,r)
    ClmnLWfluxConvCSaz = az(lonbox,latbox,ClmnLWfluxConvCS,clon,clat,r)
    ClmnSWfluxConvaz = az(lonbox,latbox,ClmnSWfluxConv,clon,clat,r)
    ClmnSWfluxConvCSaz = az(lonbox,latbox,ClmnSWfluxConvCS,clon,clat,r)
    ClmnRadfluxConvaz = az(lonbox,latbox,ClmnRadfluxConv,clon,clat,r)
    ClmnRadfluxConvCSaz = az(lonbox,latbox,ClmnRadfluxConvCS,clon,clat,r)
    SfcMoistEnthalpyFluxaz = az(lonbox,latbox,SfcMoistEnthalpyFlux,clon,clat,r)

    #Getting 10m wind speed box
    traw = int(abs(starttime-datetime(year,mo,d,hr,0,0,0)).total_seconds() / 3600)
    monthstr = time[5:7]
    WindFile = xr.open_dataset(era5MainDir+WindboxDir+str(year)+'-'+monthstr+'.nc', decode_times=False, use_cftime=False)
    #Times for file
    tarr = WindFile.time
    itarr = pd.Index(tarr)
    itlst = itarr.tolist()
    #The indexing
    tind = itlst.index(traw)
    uleft = np.array(WindFile['u10'].isel(time=tind, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
    uright = np.array(WindFile['u10'].isel(time=tind, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
    u = np.concatenate([uleft,uright],axis=1)
    vleft = np.array(WindFile['v10'].isel(time=tind, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
    vright = np.array(WindFile['v10'].isel(time=tind, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
    v = np.concatenate([vleft,vright],axis=1)
    windspeedbox = np.sqrt(u**2 + v**2)
    WindFile.close()
    windspeedboxaz = az(lonbox,latbox,windspeedbox,clon,clat,r)

    #Total Precip
    #Case of hr 6 or 7 of first of September
    if((hr==6 or hr==7) and mo==9 and d==1):
        hr1=hr-1
        d1=d
        #Get the prior month for the first file
        mo1 = mo
        mo1str = '0'+str(mo1-1)
    
        #Get the current month for the first file
        mo2 = mo
        mo2str = '0'+str(mo2)

        filetime1 = str(year)+'-'+mo1str+'.nc'
        filetime2 = str(year)+'-'+mo2str+'.nc'

        #Open the flux files
        TP_file1 = xr.open_dataset(era5MainDir+PrecipDir+filetime1, decode_times=False, use_cftime=False)
        TP_file2 = xr.open_dataset(era5MainDir+PrecipDir+filetime2, decode_times=False, use_cftime=False)
        #Get the times used for the time indexing
        t1 = int(abs(starttime-datetime(year,mo1,d1,hr1,0,0,0)).total_seconds() / 3600)
        t2 = int(abs(starttime-datetime(year,mo2,1,hr+1,0,0,0)).total_seconds() / 3600)
        #Times for file 1
        tarray1 = TP_file1.time
        itarray1 = pd.Index(tarray1)
        itlist1 = itarray1.tolist()
        #Times for file 2
        tarray2 = TP_file2.time
        itarray2 = pd.Index(tarray2)
        itlist2 = itarray2.tolist()
        #The indexing
        tind1 = itlist1.index(t1)
        tind2 = itlist2.index(t2)
        tpleft1 = np.array(TP_file1['tp'].isel(time=tind1, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
        tpright1 = np.array(TP_file1['tp'].isel(time=tind1, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
        tp1 = np.concatenate([tpleft1,tpright1],axis=1)
        tpleft2 = np.array(TP_file2['tp'].isel(time=tind2, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
        tpright2 = np.array(TP_file2['tp'].isel(time=tind2, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
        tp2 = np.concatenate([tpleft2,tpright2],axis=1)
        tp = (tp1 + tp2) / (3600*2)
        TP_file1.close()
        TP_file2.close()
        tpaz = az(lonbox,latbox,tp,clon,clat,r)
    elif((hr==0 or hr==1 or hr==2 or hr==3 or hr==4 or hr==5) and d==1):
        if(hr==0):
            hr1=23
            hr2 = hr+1
            d1=31
            mo1=mo-1
            Mostr = '0'+str(mo1)
        else:
            mo1=mo
            d1 = d
            hr1 = hr-1
            hr2 = hr+1
            Mostr = '0'+str(mo-1)
        
        filetime = str(year)+'-'+Mostr+'.nc'
        TP_file = xr.open_dataset(era5MainDir+PrecipDir+filetime, decode_times=False, use_cftime=False)
        #Get the times used for the time indexing
        t1 = int(abs(starttime-datetime(year,mo1,d1,hr1,0,0,0)).total_seconds() / 3600)
        t2 = int(abs(starttime-datetime(year,mo,d,hr2,0,0,0)).total_seconds() / 3600)
        #Times for file
        tarray = TP_file.time
        itarray = pd.Index(tarray)
        itlist = itarray.tolist()
        #The indexing
        tind1 = itlist.index(t1)
        tind2 = itlist.index(t2)
        tpleft1 = np.array(TP_file['tp'].isel(time=tind1, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
        tpright1 = np.array(TP_file['tp'].isel(time=tind1, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
        tp1 = np.concatenate([tpleft1,tpright1],axis=1)
        tpleft2 = np.array(TP_file['tp'].isel(time=tind2, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
        tpright2 = np.array(TP_file['tp'].isel(time=tind2, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
        tp2 = np.concatenate([tpleft2,tpright2],axis=1)
        tp = (tp1 + tp2) / (3600*2)
        TP_file.close()
        tpaz = az(lonbox,latbox,tp,clon,clat,r)

    else:
        if(hr==0):
            hr1 = 23
            hr2 = hr+1
            d1 = d - 1
            d2 = d
            mo1 = mo
            mo2 = mo
        elif(hr==23 and d==31):
            hr1 = hr-1
            hr2 = 1
            d2 = 1
            d1 = d
            mo1 = mo
            mo2 = mo+1
        elif(hr==23 and d!=31):
            hr1 = hr-1
            hr2 = 1
            d2 = d + 1
            d1 = d
            mo1 = mo
            mo2 = mo
        else:
            hr1 = hr-1
            hr2 = hr+1
            d1 = d
            d2 = d
            mo1 = mo
            mo2 = mo

        Mostr = '0'+str(mo)
        
        filetime = str(year)+'-'+Mostr+'.nc'
        TP_file = xr.open_dataset(era5MainDir+PrecipDir+filetime, decode_times=False, use_cftime=False)
        #Get the times used for the time indexing
        t1 = int(abs(starttime-datetime(year,mo1,d1,hr1,0,0,0)).total_seconds() / 3600)
        t2 = int(abs(starttime-datetime(year,mo2,d2,hr2,0,0,0)).total_seconds() / 3600)
        #Times for file
        tarray = TP_file.time
        itarray = pd.Index(tarray)
        itlist = itarray.tolist()
        #The indexing
        tind1 = itlist.index(t1)
        tind2 = itlist.index(t2)
        tpleft1 = np.array(TP_file['tp'].isel(time=tind1, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
        tpright1 = np.array(TP_file['tp'].isel(time=tind1, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
        tp1 = np.concatenate([tpleft1,tpright1],axis=1)
        tpleft2 = np.array(TP_file['tp'].isel(time=tind2, latitude=slice(latmin,latmax), longitude=slice(leftlonmin,leftlonmax)))
        tpright2 = np.array(TP_file['tp'].isel(time=tind2, latitude=slice(latmin,latmax), longitude=slice(rightlonmin,rightlonmax)))
        tp2 = np.concatenate([tpleft2,tpright2],axis=1)
        tp = (tp1 + tp2) / (3600*2)
        TP_file.close()
        tpaz = az(lonbox,latbox,tp,clon,clat,r)

    #Outgoing Longwave Radiation (OLR):
    OLR = NetLWToa
    OLRaz = az(lonbox,latbox,NetLWToa,clon,clat,r)

    #MSE Budget Variables Calculations
    LWavg = boxavg(ClmnLWfluxConv, latbox, lonbox)
    LWanom = ClmnLWfluxConv-LWavg

    LWcsavg = boxavg(ClmnLWfluxConvCS, latbox, lonbox)
    LWcsanom = ClmnLWfluxConvCS-LWcsavg

    OLRavg = boxavg(OLR, latbox, lonbox)
    OLRanom = OLR-OLRavg

    SWavg = boxavg(ClmnSWfluxConv, latbox, lonbox)
    SWanom = ClmnSWfluxConv-SWavg

    SWcsavg = boxavg(ClmnSWfluxConvCS, latbox, lonbox)
    SWcsanom = ClmnSWfluxConvCS-SWcsavg

    RADavg = boxavg(ClmnRadfluxConv, latbox, lonbox)
    RADanom = ClmnRadfluxConv-RADavg

    RADcsavg = boxavg(ClmnRadfluxConvCS, latbox, lonbox)
    RADcsanom = ClmnRadfluxConvCS-RADavg

    SEFavg = boxavg(SfcMoistEnthalpyFlux, latbox, lonbox)
    SEFanom = SfcMoistEnthalpyFlux-SEFavg

    HFLSavg = boxavg(hfls, latbox, lonbox)
    HFLSanom = hfls-HFLSavg

    HFSSavg = boxavg(hfss, latbox, lonbox)
    HFSSanom = hfss-HFSSavg

    #Now save the data variables to its corresponding save name created in outer loop and add the land-sea mask to convert >20% land grids to NaN
    #4D Variables
    totprecip_azsave[t,0:nslice,0:nr] = tpaz
    windbox_azsave[t,0:nslice,0:nr] = windspeedboxaz
    NetLWSfc_azsave[t,0:nslice,0:nr] = NetLWsfcaz
    NetLWSfcCS_azsave[t,0:nslice,0:nr] = NetLWsfcCSaz
    NetLWToa_azsave[t,0:nslice,0:nr] = NetLWToaaz
    NetLWToaCS_azsave[t,0:nslice,0:nr] = NetLWToaCSaz
    ClmnLWfluxConv_azsave[t,0:nslice,0:nr] = ClmnLWfluxConvaz
    ClmnLWfluxConvCS_azsave[t,0:nslice,0:nr] = ClmnLWfluxConvCSaz
    NetSWSfc_azsave[t,0:nslice,0:nr] = NetSWsfcaz
    NetSWSfcCS_azsave[t,0:nslice,0:nr] = NetSWsfcCSaz
    NetSWToa_azsave[t,0:nslice,0:nr] = NetSWToaaz
    NetSWToaCS_azsave[t,0:nslice,0:nr] = NetSWToaCSaz
    ClmnSWfluxConv_azsave[t,0:nslice,0:nr] = ClmnSWfluxConvaz
    ClmnSWfluxConvCS_azsave[t,0:nslice,0:nr] = ClmnSWfluxConvCSaz
    ClmnRadfluxConv_azsave[t,0:nslice,0:nr] = ClmnRadfluxConvaz
    ClmnRadfluxConvCS_azsave[t,0:nslice,0:nr] = ClmnRadfluxConvCSaz
    OLR_azsave[t,0:nslice,0:nr] = OLRaz
    hfls_azsave[t,0:nslice,0:nr] = hflsaz
    hfss_azsave[t,0:nslice,0:nr] = hfssaz
    sfcMoistEnthalpyFlux_azsave[t,0:nslice,0:nr] = SfcMoistEnthalpyFluxaz
    #Azmean Vars
    LHF_azmean[t,:] = azmean(lonbox,latbox,(hfls),clon,clat,r)
    SHF_azmean[t,:] = azmean(lonbox,latbox,(hfss),clon,clat,r)
    #CpTdiff_azmean[s,t,:] = azmean(data.longitude[s,t,:],data.latitude[s,t,:],data.CpTDiff[s,t,:,:],data.centerLon[s,t],data.centerLat[s,t],r)
    #LvQdiff_azmean[s,t,:] = azmean(data.longitude[s,t,:],data.latitude[s,t,:],data.LvQDiff[s,t,:,:],data.centerLon[s,t],data.centerLat[s,t],r)
    #Tdiff_azmean[s,t,:] = azmean(data.longitude[s,t,:],data.latitude[s,t,:],data.TDiff[s,t,:,:],data.centerLon[s,t],data.centerLat[s,t],r)
    #Qdiff_azmean[s,t,:] = azmean(data.longitude[s,t,:],data.latitude[s,t,:],data.QDiff[s,t,:,:],data.centerLon[s,t],data.centerLat[s,t],r)
    wind_azmean[t,:] = azmean(lonbox,latbox,(windspeedbox),clon,clat,r)
    totprecip_azmean[t,:] = azmean(lonbox,latbox,(tp),clon,clat,r)
    #2D Variables
    Clat_save[t] = clat
    Clon_save[t] = clon
    year_save[t] = int(year)
    month_save[t] = mo
    day_save[t] = d
    hour_save[t] = hr

##### Save the variables, regular variables for each year and budget variables for each year
regvars_ds = xr.Dataset(
data_vars = dict(
    TotPrecipratecirc=(['numsteps','nslice','numr'],totprecip_azsave,{'units':'m/s','long_name':'Total Precip per second Box','_FillValue':-9999,'GridType':'Radial Grid'}),
    windcirc=(['numsteps','nslice','numr'],windbox_azsave,{'units':'m/s','long_name':'10m Wind Speed Box','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetLWSfc=(['numsteps','nslice','numr'],NetLWSfc_azsave,{'units':'W/m^2','long_name':'Net LW Sfc','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetLWSfcCS=(['numsteps','nslice','numr'],NetLWSfcCS_azsave,{'units':'W/m^2','long_name':'Net LW Sfc CS','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetLWToa=(['numsteps','nslice','numr'],NetLWToa_azsave,{'units':'W/m^2','long_name':'Net LW TOA','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetLWToaCS=(['numsteps','nslice','numr'],NetLWToaCS_azsave,{'units':'W/m^2','long_name':'Net LW TOA CS','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetSWSfc=(['numsteps','nslice','numr'],NetSWSfc_azsave,{'units':'W/m^2','long_name':'Net SW Sfc','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetSWSfcCS=(['numsteps','nslice','numr'],NetSWSfcCS_azsave,{'units':'W/m^2','long_name':'Net SW Sfc CS','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetSWToa=(['numsteps','nslice','numr'],NetSWToa_azsave,{'units':'W/m^2','long_name':'Net SW TOA','_FillValue':-9999,'GridType':'Radial Grid'}),
    NetSWToaCS=(['numsteps','nslice','numr'],NetSWToaCS_azsave,{'units':'W/m^2','long_name':'Net SW TOA CS','_FillValue':-9999,'GridType':'Radial Grid'}),
    ClmnLWfluxConv=(['numsteps','nslice','numr'],ClmnLWfluxConv_azsave,{'units':'W/m^2','long_name':'Column LW Flux Convergence','_FillValue':-9999,'GridType':'Radial Grid'}),
    ClmnLWfluxConvCS=(['numsteps','nslice','numr'],ClmnLWfluxConvCS_azsave,{'units':'W/m^2','long_name':'Column LW Flux Convergence CS','_FillValue':-9999,'GridType':'Radial Grid'}),
    ClmnSWfluxConv=(['numsteps','nslice','numr'],ClmnSWfluxConv_azsave,{'units':'W/m^2','long_name':'Column SW Flux Convergence','_FillValue':-9999,'GridType':'Radial Grid'}),
    ClmnSWfluxConvCS=(['numsteps','nslice','numr'],ClmnSWfluxConvCS_azsave,{'units':'W/m^2','long_name':'Column SW Flux Convergence CS','_FillValue':-9999,'GridType':'Radial Grid'}),
    ClmnRadfluxConv=(['numsteps','nslice','numr'],ClmnRadfluxConv_azsave,{'units':'W/m^2','long_name':'Column Radiative Flux Convergence','_FillValue':-9999,'GridType':'Radial Grid'}),
    ClmnRadfluxConvCS=(['numsteps','nslice','numr'],ClmnRadfluxConvCS_azsave,{'units':'W/m^2','long_name':'Column Radiative Flux Convergence','_FillValue':-9999,'GridType':'Radial Grid'}),
    OLR=(['numsteps','nslice','numr'],OLR_azsave,{'units':'W/m^2','long_name':'Outgoing LW Radiation','_FillValue':-9999,'GridType':'Radial Grid'}),
    hfls=(['numsteps','nslice','numr'],hfls_azsave,{'units':'W/m^2','long_name':'Surface Upward Latent Heat Flux','_FillValue':-9999,'GridType':'Radial Grid'}),
    hfss=(['numsteps','nslice','numr'],hfss_azsave,{'units':'W/m^2','long_name':'Surface Upward Sensible Heat Flux','_FillValue':-9999,'GridType':'Radial Grid'}),
    SEF=(['numsteps','nslice','numr'],sfcMoistEnthalpyFlux_azsave,{'units':'W/m^2','long_name':'Surface Moist Enthalpy Flux','_FillValue':-9999,'GridType':'Radial Grid'}),
    Azmean_LHF=(['numsteps','numr'],LHF_azmean,{'units':'W*m^-2','long_name':'Az. Mean of Latent heat flux','_FillValue':-9999,'GridType':'N/A'}),
    Azmean_SHF=(['numsteps','numr'],SHF_azmean,{'units':'W*m^-2','long_name':'Az. Mean of Sensible heat flux','_FillValue':-9999,'GridType':'N/A'}),
    Azmean_Wind=(['numsteps','numr'],wind_azmean,{'units':'m/s','long_name':'Az. Mean of 10m Wind Speed','_FillValue':-9999,'GridType':'N/A'}),
    Azmean_TP=(['numsteps','numr'],totprecip_azmean,{'units':'m','long_name':'Az. Mean of Tot. Precip','_FillValue':-9999,'GridType':'N/A'}),
    wind_azmean=(['numsteps','numr'],windbin_azmean,{'units':'m/s','long_name':'Azmean dim. of Maximum Wind Speed','_FillValue':-9999,'GridType':'N/A'}),
    maxwind=(['numsteps'],maxwind_save,{'units':'m/s','long_name':'Maximum Wind Speed','_FillValue':-9999,'GridType':'Lat/Lon Grid'}),
    minSLP=(['numsteps'],minSLP_save,{'units':'hPa','long_name':'Minimum Sea Level Pressure','_FillValue':-9999,'GridType':'Lat/Lon Grid'}),
    centerLat=(['numsteps'],Clat_save,{'units':'Degrees','long_name':'TC Center Latitude Position','_FillValue':-9999,'GridType':'0.25 deg Latitude Spacing'}),
    centerLon=(['numsteps'],Clon_save,{'units':'Degrees','long_name':'TC Center Longitude Position','_FillValue':-9999,'GridType':'0.25 deg Longitude Spacing'}),
    year=(['numsteps'],year_save,{'units':'Year of given storm','long_name':'year'}),
    month=(['numsteps'],month_save,{'units':'Month of given storm','long_name':'month'}),
    day=(['numsteps'],day_save,{'units':'Day of given storm','long_name':'day'}),
    hour=(['numsteps'],hour_save,{'units':'Hour of given storm','long_name':'hour'})
)
)
regvars_ds.to_netcdf('ERA-5_CampaignTrack_Regular_Variables_Hourly_AzimuthalProfiles_NoLSmask_'+str(year)+'.nc')
regvars_ds.close()