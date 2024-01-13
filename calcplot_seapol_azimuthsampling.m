%calculate SEA-POL sampling of azimuths
clear all

azimuths = 0:5:360; %WRT N
hoursTOT = 1064; %from python

%read in coverage by azimuth calculated in python
coverage = [ 575.,  441.,  441.,  441.,  441.,  441.,  806.,  859.,  859.,        859.,  859.,  859.,  859.,  859.,  859.,  859.,  859.,  859.,        821.,  437.,  437.,  437.,  437.,  508.,  642.,  642.,  642.,        642.,  642.,  642.,  642.,  642.,  642.,  642.,  642.,  642.,        642.,  642.,  642.,  642.,  642.,  680., 1064., 1064., 1064.,       1064., 1064., 1045., 1045., 1045., 1045., 1045., 1045., 1045.,       1045.,  680.,  627.,  627.,  627.,  627.,  627.,  627.,  627.,        627.,  627.,  627.,  627.,  627.,  627.,  627.,  646.,  646.,        575.];

figure;
plot(azimuths,coverage,'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('Azimuth (degrees from N)')
ylabel('Hours')
title('SEA-POL Sampling')
xlim([0 360])
gcfsavepdf('seapolsampling-hours.pdf')

%polar plot of relative frequency
radius = 0:50:150;
coverage2 = repmat(coverage,length(radius),1);
figure;
[h1,l1] = polarPcolor(radius,azimuths,coverage2./hoursTOT,'ncolor',12,'Nspokes',9,'Ncircles',4);
set(gca,'FontSize',16)
colormap(gca,flipud(copper(12)))
caxis([0.4 1.0])
title('SEA-POL Sampling')
gcfsavepdf('seapolsampling-percent.pdf')