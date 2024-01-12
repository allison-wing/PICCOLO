%calculate SEA-POL sampling of azimuths
clear all

azimuths = 0:5:360; %WRT N
%hours of Meteor heading at each azimuth (from python script)
hoursN = 422;
hoursS = 418;
hoursW = 134; 
hoursSW = 71; 
hoursSE = 19;

%add up hours of coverage at each azimuth (WRT N) given that SEAPOL is
%blanked 90-210 degrees WRT ship heading (that is, SEA-POL samples 0-90 and
%210-360 WRT N when ship heading is N)
coverage = zeros(1,73);

%N
i1=find(azimuths==90);
i2= find(azimuths==210);
coverage(1:i1)=hoursN;
coverage(i2:end)=hoursN;

%S
i1=find(azimuths==30);
i2 = find(azimuths==270);
coverage(i1:i2) = coverage(i1:i2) + hoursS;

%W
i1=find(azimuths==120);
coverage(i1:end) = coverage(i1:end) + hoursW;

%SW
i1=find(azimuths==75);
i2=find(azimuths==315);
coverage(i1:i2) = coverage(i1:i2) +hoursSW;

%SE
i1=find(azimuths==135);
coverage(1:i1) = coverage(1:i1) + hoursSE;
i2=find(azimuths==255);
coverage(i2:end) = coverage(i2:end) + hoursSE;

figure;
plot(azimuths,coverage,'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('Azimuth (degrees from N)')
ylabel('Hours')
title('SEA-POL Sampling by Azimuth')
xlim([0 360])
% gcfsavepdf('seapolsampling-hours.pdf')

%polar plot of relative frequency
hoursTOT = hoursN + hoursS + hoursW + hoursSW + hoursSE;
radius = 0:50:150;
coverage2 = repmat(coverage,length(radius),1);
figure; 
[h1,l1] = polarPcolor(radius,azimuths,coverage2./hoursTOT,'ncolor',12,'Nspokes',9,'Ncircles',4);
set(gca,'FontSize',16)
colormap(gca,flipud(copper(12)))
caxis([0.4 1.0])
title('SEA-POL Sampling by Azimuth')
% gcfsavepdf('seapolsampling-percent.pdf')