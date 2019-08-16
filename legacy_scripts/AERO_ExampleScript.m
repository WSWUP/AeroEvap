%C. Pearson 11/09/2016
%Runs Aero Function On Example Dataset (TESTDATA.MAT)
%Test Data contains 30 minute data from Lake Lahontan Buoy

load TESTDATA.mat

%Calculate E (mm/30min); Ce; VPD (vapor pressure deficit) using aero
%function; 2 meter sensor height; 1800 sec= 30 min

for i=1:length(TESTDATA.DATENUM)
    [E(i),Ce(i),VPD(i)]=aero(TESTDATA.DATENUM(i),TESTDATA.WIND_MS(i),TESTDATA.BP_MB(i),...
        TESTDATA.AIRTEMP_C(i),TESTDATA.SKINTEMP_C(i),TESTDATA.RH(i),...
        2,1800);
end

%% Plot Data used in aero function
figure
ax1=subplot(5,1,4)
plot(TESTDATA.DATENUM,TESTDATA.BP_MB)
title('Barometric Pressure')
ylabel('mb')

ax2=subplot(5,1,5)
plot(TESTDATA.DATENUM,E)
ylabel('mm/30min')
title('30 Min Evaporation')

ax3=subplot(5,1,2)
plot(TESTDATA.DATENUM,TESTDATA.RH)
title('Relative Humidity')
ylabel('Percent')
linkaxes([ax1 ax2 ax3],'x')

ax4=subplot(5,1,1)
plot(TESTDATA.DATENUM,TESTDATA.WIND_MS)
title('Windspeed')
ylabel('m/s')


ax5=subplot(5,1,3)
plot(TESTDATA.DATENUM,TESTDATA.SKINTEMP_C,'r')
hold on
plot(TESTDATA.DATENUM,TESTDATA.AIRTEMP_C,'b')
legend('Skin Temp','Air Temp')
ylabel('C')
title('Air and Water Temp')

% ax2=subplot(5,1,5)
% plot(TESTDATA.DATENUM,TESTDATA.E)
% ylabel('mm/30min')
% title('30 Min Evaporation')


linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
dynamicDateTicks([ax1 ax2 ax3 ax4 ax5],'linked')
