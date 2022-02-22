%% Demonstrate how to use my processing code
close all
clear 

cd ~/Documents/data/OpTrap

% Path for data:
dPath = '2022_02_01/bead_exp02';

% Load data and images
data = bead_loadData(dPath, true);

% Pre-process data
data = bead_preProcessCentres(data);

% Plot position-time data in physical units
bead_plotProData(data);
% Show data on top of full field-of-view image
bead_imageAndScatter(data);

% Calculate MSD
data = bead_normMSD(data);

% Show the FFT of position data
data = bead_fft_scaled(data, true);

% Filter out noise peaks
data.opts.bandstop = [92 93; 124 125; 152 153];
data = bead_filter_bandstop(data);

% Inspect filtered data
data = bead_fft_scaled(data,true);

% Recalculate MSD
data = bead_normMSD(data, 'forceRun', true);

% Quantify MSD
tRanges = repmat({{[1e-4 8e-3] [0.1 5]}},1,2);
[cTau, fps, RMSE] = msd_cornerator(data.pro.amsdObj, 0, tRanges);
