%% Anna and Megan 

% Instructions: Follow through this code step by step, while also referring
% to the overall instructions and questions from the lab assignment sheet.

%% 1. Read in the monthly gridded CO2 data from the .csv file
% The data file is included in your repository as �LDEO_GriddedCO2_month_flux_2006c.csv�
% Your task is to write code to read this in to MATLAB
% Hint: you can again use the function �readtable�, and use your first data lab code as an example.
filename = 'LDEO_GriddedCO2_month_flux_2006c.csv'
CO2data = readtable(filename);

%% 2a. Create new 3-dimensional arrays to hold reshaped data
%Find each unique longitude, latitude, and month value that will define
%your 3-dimensional grid
longrid = unique(CO2data.LON); %finds all unique longitude values
latgrid = unique(CO2data.LAT);%<-- following the same approach, find all unique latitude values
monthgrid = unique(CO2data.MONTH); %<-- following the same approach, find all unique months

%Create empty 3-dimensional arrays of NaN values to hold your reshaped data
    %You can make these for any variables you want to extract - for this
    %lab you will need PCO2_SW (seawater pCO2) and SST (sea surface
    %temperature)
SWco2 = NaN*zeros(length(longrid), length(latgrid), length(monthgrid));
SStemp = NaN*zeros(length(longrid), length(latgrid), length(monthgrid));


%% 2b. Pull out the seawater pCO2 (PCO2_SW) and sea surface temperature (SST)
%data and reshape it into your new 3-dimensional arrays

for i = 1:length(CO2data.LON)
    
   lonind = find(CO2data.LON(i) == longrid);
   latind = find(CO2data.LAT(i) == latgrid);
   monthind = find(CO2data.MONTH(i) == monthgrid);
          
   SWco2(lonind, latind, monthind) = CO2data.PCO2_SW(i);
   
   % adding to SST
          
   SStemp(lonind, latind, monthind) = CO2data.SST(i);
   
   % for part 5- adding in atmospheric pCO2
   
   AIRco2(lonind,latind, monthind) = CO2data.PCO2_AIR(i);

end

%% 3a. Make a quick plot to check that your reshaped data looks reasonable
%Use the imagesc plotting function, which will show a different color for
%each grid cell in your map. Since you can't plot all months at once, you
%will have to pick one at a time to check - i.e. this example is just for
%January

%January
%imagesc(SStemp(:,:,1))

%July
imagesc(SStemp(:,:,7))
%% 3b. Now pretty global maps of one month of each of SST and pCO2 data.
%I have provided example code for plotting January sea surface temperature
%(though you may need to make modifications based on differences in how you
%set up or named your variables above).

figure(1); clf
worldmap world
contourfm(latgrid, longrid, SStemp(:,:,1)','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('January Sea Surface Temperature (^oC)')
%%
%Check that you can make a similar type of global map for another month
%and/or for pCO2 using this approach. Check the documentation and see
%whether you can modify features of this map such as the contouring
%interval, color of the contour lines, labels, etc.

figure(2); clf
worldmap world
contourfm(latgrid, longrid, SStemp(:,:,7)','linecolor','black');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('July Sea Surface Temperature (^oC)')

%% 4. Calculate and plot a global map of annual mean pCO2

%pCO2mean = mean(CO2data.PCO2_SW);

% annual mean seawater pCO2 for each location 
aPCO2 = mean(SWco2,3)';

%mapping annual mean pCO2

figure(3); clf
worldmap world
contourfm(latgrid, longrid, aPCO2,'linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Annual Mean pCO2 (^oC)')

%% 5. Calculate and plot a global map of the difference between the annual mean seawater and atmosphere pCO2
% used atmospheric pCO2 data from NOAA website
%addpath(C:\Users\megan\Documents\GitHub\Common-functions)
monthmeanair = mean(AIRco2, 3); % annual mean Atmos. pCO2 for each location

figure(4); clf
worldmap world
contourfm(latgrid, longrid, aPCO2-368.84,'linecolor','none');
geoshow('landareas.shp','FaceColor','black')
title('Annual Mean Difference Between Seawater and Atmospheric pCO2 (^oC)')
cmocean('balance', 'pivot',0);
colorbar

%% 6. Calculate relative roles of temperature and of biology/physics in controlling seasonal cycle
%variables for equation 1
aTemp = mean(SStemp,3); %annual mean temp for each location- now a 2D array
aPCO2 = mean(SWco2,3); %annual mean seawater pCO2 for each location- now a 2D array

%equation 1 (biophysical effects)
aTemp1 = repmat(aTemp,1,1,12); %repmat to change aTemp from 2D to 3D
eq1 = SWco2.*exp(0.0423.*(aTemp1-SStemp)); 
pCO2_BP = eq1; %3D array of BP effects

%equation 2 (temp effects)
eq2 = aPCO2.*exp(0.0423.*(SStemp-aTemp1)); %aTemp1 is in here now
pCO2_T = eq2; %3D array of temp effects

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)

% figuring out lat&lon

latstn = latgrid-4 %ends at 76 [-80, -76, -72...72, 76]
lonstn = longrid - 2.5 %ends at 355 [0, 5, 10, 15...355] 
% nevermind we don't need this (above)

%[minval,minID]=min(A, B) = A is min, B is where it is
%diff between latgrid and the actual lat of station

%BATS lat lon
[minval1, id1] = min(abs(latgrid - 32)) %gives nearest LAT integer & where it is
[minval11, id11] = min(abs(longrid - 296)) %gives nearest LON integer & where it is

%Papa lat lon
[minval2, id2] = min(abs(latgrid - 50))
[minval22, id22] = min(abs(longrid - 215))

%ROSS lat lon
[minval3, id3] = min(abs(latgrid - -76))
[minval33, id33] = min(abs(longrid - 183))

%% STATION PLOTS FOR PART 7
%%%%%%%%%%%%
%Ocean station BATS 32^0 50' N, -64^0 10' W 

figure(5); clf
plot(squeeze(SStemp(id11,id1,:))) %sea surface temp (id11,id1,:)

hold on
plot(squeeze(SWco2(id11,id1,:))) %seawater pCO2

hold on
plot(squeeze(pCO2_T(id11,id1,:))) %seawater pCO2 temp effects

hold on
plot(squeeze(pCO2_BP(id11,id1,:))) %seawater pCO2 temp effects

xlim([0 12])
ylim([0 500])
legend('SST', 'CO2 Levels', 'Temp Effect','BP Effect', 'location', 'east')
title('Seasonal Cycle For BATS')
xlabel('Months')
ylabel('pCO2')
hold off



%%
%%%%%%%%%%%%%%%%
%Ocean Station Papa (50 N lat / -145 W lon)
figure(6); clf
plot(squeeze(SStemp(id22,id2,:)))

hold on
plot(squeeze(SWco2(id22, id2,:)))
% 
hold on
plot(squeeze(pCO2_T(id22,id2,:)))
% 
hold on
plot(squeeze(pCO2_BP(id22,id2,:)))

legend('SST', 'CO2 Levels', 'Temp Effect','BP Effect', 'location', 'east')
title('Seasonal Cycle For Papa')
xlabel('Months')
ylabel('pCO2')
hold off
%% 
%%%%%%%%%%%%%%%%%%ROSS is the only one that works for me and I don't know why!
%ROSS (76 S lat / -177.5 W lon)
figure(7); clf
plot(squeeze(SStemp(id33,id3,:)))

hold on
plot(squeeze(SWco2(id33,id3,:)))
% 
hold on
plot(squeeze(pCO2_T(id33,id3,:)))
% 
hold on
plot(squeeze(pCO2_BP(id33,id3,:)))

legend('SST', 'CO2 Levels', 'Temp Effect','BP Effect', 'location', 'east')
title('Seasonal Cycle For ROSS')
xlabel('Months')
ylabel('pCO2')
hold off
%End of stations

%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on these maps the locations of the three stations for which you plotted the
% seasonal cycle above


%%%%% note I had to use the ' on each effect in global map codes

%page 1614 in Takahashi paper 

%First map is biophysical pCO2 effect seasonal amplitude 
%equation in paper : (change in pco2)bio = (pco2 at Tmean)max - (pco2 at
%Tmean)min

%EQUATION CODE 
cpco2Bio = max(pCO2_BP, [], 3) - min(pCO2_BP, [], 3); %now is 2D

%GLOBAL MAP CODE
figure(8); clf
worldmap world
contourfm(latgrid, longrid, cpco2Bio','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Biophysical Effects on pCO2 Seasonal Amplitude')
plotm(latgrid(id1), longrid(id11), 'Marker', '.', 'Color', 'yellow', 'MarkerSize', 40);
plotm(latgrid(id2), longrid(id22),'Marker', '.', 'Color', 'red', 'MarkerSize', 40);
plotm(latgrid(id3), longrid(id33),'Marker', '.', 'Color', 'magenta', 'MarkerSize', 40);
legend('','BATS', 'Papa', 'Ross', 'Location', 'SouthEastOutside')

%%

%Second map is temperature amplitude
%equation in paper: (change in pco2)temp = (pco2 at Tobs)max - (pco2 at
%Tobs) min

%EQUATION CODE
cpco2Temp = max(pCO2_T, [], 3) - min(pCO2_T, [], 3);

%GLOBAL MAP CODE
figure(9); clf
worldmap world
contourfm(latgrid, longrid, cpco2Temp','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Temperature Effect Seasonal Amplitude')
plotm(latgrid(id1), longrid(id11), 'Marker', '.', 'Color', 'yellow', 'MarkerSize', 40);
plotm(latgrid(id2), longrid(id22),'Marker', '.', 'Color', 'red', 'MarkerSize', 40);
plotm(latgrid(id3), longrid(id33),'Marker', '.', 'Color', 'magenta', 'MarkerSize', 40);
legend('','BATS', 'Papa', 'Ross', 'Location', 'SouthEastOutside')

%%

%Third map is the difference between the two (T-B)
%equation in paper: (T-B) = (change in pco2)temp - (change in pco2)bio

%EQUATION CODE
diff = cpco2Temp - cpco2Bio;

%GLOBAL MAP CODE
figure(10); clf
worldmap world
contourfm(latgrid, longrid, diff','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('The Difference Between Temperature & Biophysical Amplitude')
cmocean('balance', 'pivot',0);
plotm(latgrid(id1), longrid(id11), 'Marker', '.', 'Color', 'yellow', 'MarkerSize', 40);
plotm(latgrid(id2), longrid(id22),'Marker', '.', 'Color', 'red', 'MarkerSize', 40);
plotm(latgrid(id3), longrid(id33),'Marker', '.', 'Color', 'magenta', 'MarkerSize', 40);
legend('','BATS', 'Papa', 'Ross', 'Location', 'SouthEastOutside')
