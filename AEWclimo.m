%Count how many AEW tracks pass through 1S-15N, 50W-23W between August 13
%and September 30, based on Quinton Lawton's ERA-5 AEW climatology
clear all

years = 1979:2021;
lat1 = -1;
lat2 = 20;
lon1 = -60;
lon2 = -20;

gridsize = 1; %grid size for track density

%ship track
shipLON=[-24.98, -23, -23, -23.4, -23.4, -23.8, -32, -32, -32.4, -32.4, -32.8, -47, -47, -47.4, -47.4, -47.4, -59.42]
shipLAT = [16.88, 14.5, -0.5, 14.5, 2.5, 7.5, 7.5, 14.5, 2.5, 14.5, 7.5, 7.5, 14.5, 5, 7.9, 14.9, 13.15]

%read in lat/lon arrays
ncid = netcdf.open('/Users/allisonwing/Dropbox/AEWtracks/AEW_data_postprocessed__0924_UPDATE_EXTENDED_1986_B1-6hr.nc');
lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latitude'));
lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'longitude'));
netcdf.close(ncid)

AEWcount = zeros(1,length(years)); %number of AEWs in our region/time window in each year
AEWlon_subset = nan(612,81,length(years));
AEWlat_subset = nan(612,81,length(years));
figure;
for i = 1:length(years)
    %Open netcdf file for this year
    fname = ['/Users/allisonwing/Dropbox/AEWtracks/AEW_data_postprocessed__0924_UPDATE_EXTENDED_' num2str(years(i)) '_B1-6hr.nc']
    ncid = netcdf.open(fname);
    AEW_lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lon')); %time x system
    AEW_lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lat')); %time x system
    time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %hours since 1900-01-01 00:00:00.0
    sz = size(AEW_lat);
    ntime(i) = sz(1);
    nsystem(i) = sz(2);
    
    %extract indices for times that fall within our window
    tstart = datenum(1900,1,1); %days since 0-Jan-0000
    ymd = datevec(time/24+tstart); %year month day for AEWtracks (%time x 6
    month = ymd(:,2);
    day = ymd(:,3);
    itimecheck = find((month==8 & day>=10) | (month==9 & day<=23)); %indexing in time (august 13 - september 30)
    
    %extract AEWdata over the times that fall within our window
    AEWlon = AEW_lon(itimecheck,:);
    AEWlat = AEW_lat(itimecheck,:);
    AEWlon_subset(itimecheck,1:nsystem(i),i) = AEWlon;
    AEWlat_subset(itimecheck,1:nsystem(i),i) = AEWlat;
    
    %loop over each system and check if it is in our region
    for j = 1:nsystem(i)
        %extract indices for lat and lon that fall within our region
        iregioncheck = find(AEWlat(:,j)>=lat1 & AEWlat(:,j)<=lat2 & AEWlon(:,j)>=lon1 & AEWlon(:,j)<=lon2);
        
        %if anything was found for that storm, add to count for this year
        if ~isempty(iregioncheck)
            AEWcount(i) = AEWcount(i) + 1;
        end
    end
    
    %plot tracks over all years year
    hold on
    plot(AEWlon_subset(:,:,i),AEWlat_subset(:,:,i),'k')
end

ylim([-1 20])
xlim([-60 -20])
set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
title('AEWs: 10 Aug - 23 Sep, 1979:2021')
gcfsavepdf('AEWtracks_all.pdf')

%plot for one year
whichyear = 2021;
iyear = find(years==whichyear);
figure;
plot(AEWlon_subset(:,:,iyear),AEWlat_subset(:,:,iyear),'k')
ylim([-1 20])
xlim([-60 -20])
set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
title(['AEWs: 10 Aug - 23 Sep ' num2str(whichyear)])
gcfsavepdf('AEWtracks_oneyear.pdf')

%Climatological average AEWcount for that time period/region (average
%AEWcount over all years)
climoAEWcount = mean(AEWcount)

%calculate track density
Xp = lon1:gridsize:lon2;
Yp = lat2:-gridsize:lat1; %reverse order to match AEW lat data
lx = length(Xp);
ly = length(Yp);

trden=zeros(ly,lx,length(years));

%loop over all years, all storms, all points along tracks [already
%restricted to storms passing through our box and our time period]
for i = 1:612 %track points
    for j = 1:81 %storms
        for k = 1:length(years) %years
            %if the storm can be assigned to 2 grid boxes, choose the
            %nearest one
            ix=find(Xp>=AEWlon_subset(i,j,k)-gridsize/2 & Xp<AEWlon_subset(i,j,k)+gridsize/2);
            if (length(ix)>1)
                diff = Xp(ix) - AEWlon_subset(i,j,k);
                [xx,id] = min(abs(diff));
                ixv = ix(id);
                clear ix;
                ix = ixv;
            end
            iy=find(Yp>=AEWlat_subset(i,j,k)-gridsize/2 & Yp<AEWlat_subset(i,j,k)+gridsize/2);
            if (length(iy)>1)
                diff = Yp(iy) - AEWlat_subset(i,j,k);
                [yy,id] = min(abs(diff));
                iyv = iy(id);
                clear iy;
                iy = iyv;
            end
            trden(iy,ix,k) = trden(iy,ix)+0.25; %because tracks given four times per day
            
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

%plot for one year
whichyear = 2021;
iyear = find(years==whichyear);
hold on
plot(AEWlon_subset(:,:,iyear),AEWlat_subset(:,:,iyear),'w')

%add ship track
plot(shipLON,shipLAT,'r','LineWidth',2)

axis equal
ylim([-1 20])
xlim([-60 -20])
set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
title('AEW Track Density 10 Aug - 23 Sep')

gcfsavepdf('AEWtracks_climo_2021.pdf')
