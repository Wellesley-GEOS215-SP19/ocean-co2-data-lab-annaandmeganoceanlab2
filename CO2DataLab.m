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

end

%% 3a. Make a quick plot to check that your reshaped data looks reasonable
%Use the imagesc plotting function, which will show a different color for
%each grid cell in your map. Since you can't plot all months at once, you
%will have to pick one at a time to check - i.e. this example is just for
%January

imagesc(SStemp(:,:,1))

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

pCO2mean = mean(CO2data.PCO2_SW)

figure(1); clf
worldmap world
contourfm(latgrid, longrid, pCO2mean(:,:,:),'linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Annual Mean pCO2 (^oC)')

%% 5. Calculate and plot a global map of the difference between the annual mean seawater and atmosphere pCO2
%<--

%% 6. Calculate relative roles of temperature and of biology/physics in controlling seasonal cycle
%<--

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)

%<--

%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above

%<--
