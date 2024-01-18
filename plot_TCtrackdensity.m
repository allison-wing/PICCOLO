%% Tropical cyclone track density plot from IBTrACS
clear all

%% set things
% restrict years
year1 = 1980;
year2 = 2023;
years = year1:year2;

%time window of interest
month1 = 8;
day1 = 10;
month2 = 9;
day2 = 24;

%region of interest
lat1 = -1;
lat2 = 20;
lon1 = -60;
lon2 = -20;

gridsize = 1; %grid size for track density
tfactor = 0.125; %0.25 for 6-hour track points, 0.125 (1/8th of day) for 3-hour track points

%ship track
shipLON=[-24.98, -23, -23, -23.4, -23.4, -23.8, -32, -32, -32.4, -32.4, -32.8, -47, -47, -47.4, -47.4, -47.4, -59.42];
shipLAT = [16.88, 14.5, -0.5, 14.5, 2.5, 7.5, 7.5, 14.5, 2.5, 14.5, 7.5, 7.5, 14.5, 5, 7.9, 14.9, 13.15];


%% read in data
ncid = netcdf.open('/Volumes/home/awing/ibtracs/IBTrACS.NA.v04r00.nc');

numobs = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'numobs')); %number of obs per system
time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %days since 1858-11-17 00:00:00
lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'usa_lat')); %USA lat position
lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'usa_lon')); %USA lat position
wind = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'usa_wind')); %max sustained wind (kts)
status = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'usa_status')); %what storm type
sscat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'usa_sshs')); %category on SSHS
%-5 = unknown, -4 = post-tropical/ET; -3 = disturbance; -2 = subtropical;
%-1 = depression; 0 = tropical storm; 1 = cat 1; 2 = cat 2; 3 = cat3;
%4 = at 4; 5 = cat5;

%sort out date
time(time<-9000000)=NaN;
tstart = datenum(1858,11,17); %days since 0-Jan-00

%pre-allocate
sz = size(lat);
nyears = length(years);
lat_subset = nan(sz(1),40,nyears);
lon_subset = nan(sz(1),40,nyears);
wind_subset = nan(sz(1),40,nyears);
sscat_subset = nan(sz(1),40,nyears);
time_subset = nan(sz(1),40,nyears);

%% Loop through storms and grab those that are in our year range, date range.
%and restrict to 0,6,12,18Z times
%Reorganize to be time x nstorm x year
stormcount = 0;
for nstorm = 1:length(numobs)
    %pull time data for that storm
    ymd = datevec(time(:,nstorm)+tstart);
    year(:,nstorm) = ymd(:,1);
    month(:,nstorm) = ymd(:,2);
    day(:,nstorm) = ymd(:,3);
    hour(:,nstorm) = ymd(:,4);
    itimecheck = find(year(:,nstorm)>=year1 & year(:,nstorm)<=year2 & ...
        (hour(:,nstorm)==0 | hour(:,nstorm)==3 | hour(:,nstorm)==6 | hour(:,nstorm)==9 | hour(:,nstorm)==12 | hour(:,nstorm)==15 | hour(:,nstorm)==18 | hour(:,nstorm)==21) & ...
        ((month(:,nstorm)==month2 & day(:,nstorm)>=day2) | (month(:,nstorm)==month2 & day(:,nstorm)<=day2))); %indexing in time
    
    %identify which storm number of the year this is
    iyear = find(years==year(1,nstorm)); %this storm
    if nstorm>1
        iyear_last = find(years==year(1,nstorm-1)); %last storm
        if years(iyear)~=years(iyear_last) %if in different year than last storm
            stormcount = 0; %reset storm count to zero
        end
    end
   
    
    % if anything was found for that storm, mark it
    if ~isempty(itimecheck)
        if ~exist('whichTC','var')
            whichTC = nstorm;
        else
            whichTC = [whichTC nstorm]; %build list of indices corresponding to storms that are in our time window
        end
        
        stormcount = stormcount+ 1;
        
        lat_subset(itimecheck,stormcount,iyear) = lat(itimecheck,nstorm);
        lon_subset(itimecheck,stormcount,iyear) = lon(itimecheck,nstorm);
        wind_subset(itimecheck,stormcount,iyear) = wind(itimecheck,nstorm);
        sscat_subset(itimecheck,stormcount,iyear) = sscat(itimecheck,nstorm);
        time_subset(itimecheck,stormcount,iyear) = time(itimecheck,nstorm);
    end
end

%plot for one year
whichyear = 2021;
iyear = find(years==whichyear);
figure;
plot(lon_subset(:,:,iyear),lat_subset(:,:,iyear),'k')
ylim([lat1 lat2])
xlim([lon1 lon2])
set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
title(['TCs: 10 Aug - 24 Sep ' num2str(whichyear)])
gcfsavepdf('TCtracks_oneyear.pdf')

%plot for all years
figure;
for iyear = 1:nyears
    hold on
    plot(lon_subset(:,:,iyear),lat_subset(:,:,iyear),'k')
end
ylim([lat1 lat2])
xlim([lon1 lon2])
set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
title('TCs: 10 Aug - 24 Sep, 1980:2023')
gcfsavepdf('TCtracks_all.pdf')

%% calculate track density
Xp = lon1:gridsize:lon2;
Yp = lat2:-gridsize:lat1; %reverse order to match TC lat data
lx = length(Xp);
ly = length(Yp);

trden=zeros(ly,lx,length(years));

%loop over all years, all storms, all points along tracks [already
%restricted to storms passing in our time period]
for i = 1:sz(1) %track points
    for j = 1:40 %storms
        for k = 1:nyears %years
            %if the storm can be assigned to 2 grid boxes, choose the
            %nearest one
            ix=find(Xp>=lon_subset(i,j,k)-gridsize/2 & Xp<lon_subset(i,j,k)+gridsize/2);
            if (length(ix)>1)
                diff = Xp(ix) - lon_subset(i,j,k);
                [xx,id] = min(abs(diff));
                ixv = ix(id);
                clear ix;
                ix = ixv;
            end
            iy=find(Yp>=lat_subset(i,j,k)-gridsize/2 & Yp<lat_subset(i,j,k)+gridsize/2);
            if (length(iy)>1)
                diff = Yp(iy) - lat_subset(i,j,k);
                [yy,id] = min(abs(diff));
                iyv = iy(id);
                clear iy;
                iy = iyv;
            end
            trden(iy,ix,k) = trden(iy,ix)+tfactor; %because tracks given four times per day
            
        end
    end
end

%climatological track density
trdencl = squeeze(mean(trden,3));
trdencl(trdencl==0)=NaN; %make NaN where zero

%plot
figure;
h=imagesc(Xp,Yp,trdencl);
set(gca,'Ydir','normal')
colormap(turbo)
set(h, 'AlphaData', ~isnan(trdencl)) %set to white where NaN (zero)
colorbar

%add plot for one year
whichyear = 2021;
iyear = find(years==whichyear);
hold on
plot(lon_subset(:,:,iyear),lat_subset(:,:,iyear),'w')

%add ship track
plot(shipLON,shipLAT,'r','LineWidth',2)

axis equal
ylim([lat1 lat2])
xlim([lon1 lon2])
set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
title('TC Track Density 10 Aug - 24 Sep')

gcfsavepdf('TCtracks_climo_2021.pdf')



