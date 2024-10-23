%Plot AEW tracks during BOWTIE cruise
%Plot climatological AEW tracks, actual tracks, and Meteor ship track 
%between August 16 and September 24, through 4S-20N,20W to 60W

clear all

years = 1979:2023;
lat1 = 4;
lat2 = 20;
lon1 = -60;
lon2 = -20;

gridsize = 1; %grid size for track density


%% Open ship track
ncid = netcdf.open('/Users/awing/Dropbox/ORCESTRA/data/meteor_meteo_dship_20240923.nc'); %minutes since 2024-08-14 00:00
shipLON = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
shipLAT = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
shipTIME = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
netcdf.close(ncid)

%% AEW climatology


%read in lat/lon arrays
% ncid = netcdf.open('/Users/awing/Dropbox/AEWtracks/AEW_data_postprocessed__0924_UPDATE_EXTENDED_1986_B1-6hr.nc');
ncid = netcdf.open('/Users/awing/Dropbox/AEWtracks/ERA5_AEW_TRACKS/AEW_tracks_post_processed_year_1986.nc');
lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latitude'));
lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'longitude'));
netcdf.close(ncid)

AEWcount = zeros(1,length(years)); %number of AEWs in our region/time window in each year
AEWlon_subset = nan(612,81,length(years));
AEWlat_subset = nan(612,81,length(years));

for i = 1:length(years)
    %Open netcdf file for this year
%     fname = ['/Users/awing/Dropbox/AEWtracks/AEW_data_postprocessed__0924_UPDATE_EXTENDED_' num2str(years(i)) '_B1-6hr.nc']
    fname = ['/Users/awing/Dropbox/AEWtracks/ERA5_AEW_TRACKS/AEW_tracks_post_processed_year_' num2str(years(i)) '.nc']
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
    itimecheck = find((month==8 & day>=16) | (month==9 & day<=24)); %indexing in time (august 16 - september 24)
    
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
end

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


%% AEW tracks

%pre-allocate for AEW tracks from GFS analysis
AEWlon_cat = nan(500,65);
AEWlat_cat = nan(500,65);
AEWtime_cat = nan(500,1);

%Open AEW tracks
tstart = datenum(1900,1,1); %days since 0-Jan-0000
%File 1
AEWfile_mon = 8;
AEWfile_day = 23;
AEWfile_hour = 6;
ncid = netcdf.open('/Users/awing/Dropbox/ORCESTRA/data/AEW/tracks_2024082306.nc');
AEWlon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lon_smooth')); %time x system
AEWlat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lat_smooth')); %time x system
AEWtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %hours since 01/01/1900

%convert times
analtime = 24*(datenum(2024,AEWfile_mon,AEWfile_day,AEWfile_hour,0,0) - tstart); %hours

%extract just the times up to the analysis time
ianaltime = find(AEWtime==analtime);
sz = size(AEWlon);
AEWlon_cat(1:ianaltime,1:sz(2)) = AEWlon(1:ianaltime,:);
AEWlat_cat(1:ianaltime,1:sz(2)) = AEWlat(1:ianaltime,:);
AEWtime_cat(1:ianaltime) = AEWtime(1:ianaltime);

%File 2
AEWfile_mon = 8;
AEWfile_day = 30;
AEWfile_hour = 6;
ncid = netcdf.open('/Users/awing/Dropbox/ORCESTRA/data/AEW/tracks_2024083006.nc');
AEWlon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lon_smooth')); %time x system
AEWlat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lat_smooth')); %time x system
AEWtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %hours since 01/01/1900

%convert times
analtime2 = 24*(datenum(2024,AEWfile_mon,AEWfile_day,AEWfile_hour,0,0) - tstart); %hours

%extract just the times up to the analysis time
ianaltime2 = find(AEWtime==analtime2);
sz2 = size(AEWlon);
AEWlon_cat(ianaltime+1:ianaltime+ianaltime2-1,sz(2)+1:sz(2)+sz2(2)) = AEWlon(2:ianaltime2,:);
AEWlat_cat(ianaltime+1:ianaltime+ianaltime2-1,sz(2)+1:sz(2)+sz2(2)) = AEWlat(2:ianaltime2,:);
AEWtime_cat(ianaltime+1:ianaltime+ianaltime2-1) = AEWtime(2:ianaltime2);

%File 3
AEWfile_mon = 9;
AEWfile_day = 6;
AEWfile_hour = 6;
ncid = netcdf.open('/Users/awing/Dropbox/ORCESTRA/data/AEW/tracks_2024090606.nc');
AEWlon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lon_smooth')); %time x system
AEWlat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lat_smooth')); %time x system
AEWtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %hours since 01/01/1900

%convert times
analtime3 = 24*(datenum(2024,AEWfile_mon,AEWfile_day,AEWfile_hour,0,0) - tstart); %hours

%extract just the times up to the analysis time
ianaltime3 = find(AEWtime==analtime3);
sz3 = size(AEWlon);
AEWlon_cat(ianaltime+ianaltime2:ianaltime+ianaltime2+ianaltime3-2,sz(2)+sz2(2):sz(2)+sz2(2)+sz3(2)-1) = AEWlon(2:ianaltime3,:);
AEWlat_cat(ianaltime+ianaltime2:ianaltime+ianaltime2+ianaltime3-2,sz(2)+sz2(2):sz(2)+sz2(2)+sz3(2)-1) = AEWlat(2:ianaltime3,:);
AEWtime_cat(ianaltime+ianaltime2:ianaltime+ianaltime2+ianaltime3-2) = AEWtime(2:ianaltime3);

%File 4
AEWfile_mon = 9;
AEWfile_day = 13;
AEWfile_hour = 12;
ncid = netcdf.open('/Users/awing/Dropbox/ORCESTRA/data/AEW/tracks_2024091312.nc');
AEWlon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lon_smooth')); %time x system
AEWlat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lat_smooth')); %time x system
AEWtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %hours since 01/01/1900

%convert times
analtime4 = 24*(datenum(2024,AEWfile_mon,AEWfile_day,AEWfile_hour,0,0) - tstart); %hours

%extract just the times up to the analysis time
ianaltime4 = find(AEWtime==analtime4);
sz4 = size(AEWlon);
AEWlon_cat(ianaltime+ianaltime2+ianaltime3-1:ianaltime+ianaltime2+ianaltime3+ianaltime4-2,sz(2)+sz2(2)+sz3(2):sz(2)+sz2(2)+sz3(2)+sz4(2)-1) = AEWlon(1:ianaltime4,:);
AEWlat_cat(ianaltime+ianaltime2+ianaltime3-1:ianaltime+ianaltime2+ianaltime3+ianaltime4-2,sz(2)+sz2(2)+sz3(2):sz(2)+sz2(2)+sz3(2)+sz4(2)-1) = AEWlat(1:ianaltime4,:);
AEWtime_cat(ianaltime+ianaltime2+ianaltime3-1:ianaltime+ianaltime2+ianaltime3+ianaltime4-2) = AEWtime(1:ianaltime4);

%File 5
AEWfile_mon = 9;
AEWfile_day = 20;
AEWfile_hour = 12;
ncid = netcdf.open('/Users/awing/Dropbox/ORCESTRA/data/AEW/tracks_2024092012.nc');
AEWlon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lon_smooth')); %time x system
AEWlat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lat_smooth')); %time x system
AEWtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %hours since 01/01/1900

%convert times
analtime5 = 24*(datenum(2024,AEWfile_mon,AEWfile_day,AEWfile_hour,0,0) - tstart); %hours

%extract just the times up to the analysis time
ianaltime5 = find(AEWtime==analtime5);
sz5 = size(AEWlon);
AEWlon_cat(ianaltime+ianaltime2+ianaltime3+ianaltime4-1:ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5-3,sz(2)+sz2(2)+sz3(2)+sz4(2):sz(2)+sz2(2)+sz3(2)+sz4(2)+sz5(2)-1) = AEWlon(2:ianaltime5,:);
AEWlat_cat(ianaltime+ianaltime2+ianaltime3+ianaltime4-1:ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5-3,sz(2)+sz2(2)+sz3(2)+sz4(2):sz(2)+sz2(2)+sz3(2)+sz4(2)+sz5(2)-1) = AEWlat(2:ianaltime5,:);
AEWtime_cat(ianaltime+ianaltime2+ianaltime3+ianaltime4-1:ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5-3) = AEWtime(2:ianaltime5);

%File 6
AEWfile_mon = 9;
AEWfile_day = 24;
AEWfile_hour = 12;
ncid = netcdf.open('/Users/awing/Dropbox/ORCESTRA/data/AEW/tracks_2024092412.nc');
AEWlon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lon_smooth')); %time x system
AEWlat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'AEW_lat_smooth')); %time x system
AEWtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')); %hours since 01/01/1900

%convert times
analtime6 = 24*(datenum(2024,AEWfile_mon,AEWfile_day,AEWfile_hour,0,0) - tstart); %hours
endttime = 24*(datenum(2024,9,24,3,0,0) - tstart);

%extract just the times up to end of campaign
ianaltime6 = find(AEWtime==endttime);
sz6 = size(AEWlon);
AEWlon_cat(ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5-1:ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5+ianaltime6-3,sz(2)+sz2(2)+sz3(2)+sz4(2):sz(2)+sz2(2)+sz3(2)+sz4(2)+sz5(2)+sz6(2)-1) = AEWlon(2:ianaltime6,:);
AEWlat_cat(ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5-1:ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5+ianaltime6-3,sz(2)+sz2(2)+sz3(2)+sz4(2):sz(2)+sz2(2)+sz3(2)+sz4(2)+sz5(2)+sz6(2)-1) = AEWlat(2:ianaltime6,:);
AEWtime_cat(ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5-1:ianaltime+ianaltime2+ianaltime3+ianaltime4+ianaltime5+ianaltime6-3) = AEWtime(2:ianaltime6);



%plot tracks
szAEW = size(AEWlon_cat);
nAEW = szAEW(2);

ymd = datevec(AEWtime_cat/24+tstart);
year = ymd(:,1);
month = ymd(:,2);
day = ymd(:,3);
hour = ymd(:,4);


for i = 1:nAEW
    hold on
    plot(AEWlon_cat(:,i),AEWlat_cat(:,i),'LineWidth',2)
end
plot(shipLON,shipLAT,'r','LineWidth',3) %whole ship track
% 
axis equal
ylim([lat1 lat2])
xlim([lon1 lon2])
set(gca,'FontSize',16)
xlabel('Longitude')
ylabel('Latitude')
title(['AEWs: 2024.08.16:06UTC - 2024.09.24:03UTC'])

gcfsavepdf('./figures/AEWtracks_wclimo.pdf')
% gcfsavepdf('./figures/AEWtracks.pdf')
