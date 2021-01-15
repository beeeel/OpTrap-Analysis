%% High-pass filter position and then calculate allan variance. 
% Do this for several frequencies and plot all the variances together

%% Setup - first run bead_processing_v0 to get relavant variables
fpasses = [0.5 1 5 10 50 100]; % Pass frequency
fs = length(xCentres)./ (1e-3 * (max(timeVec) - min(timeVec))); % Sampling frequency in Hz
cropT = cropTs{fileIdx}; % Frame range to take (decided by eying the data)
cropTHPval = 1; % Number of frames to crop to remove ringing caused by HP

%% Preallocate 
% Crop indices for after high-passing (remove ringing)
cropTHP = [cropTHPval, diff(cropT) + 1 - cropTHPval];
% Number of points in final HP vector
nPoints = diff(cropT) - 2 * cropTHPval + 2;
% Lag times
lagTimes = calcLagTimes(nPoints);
YHP = zeros(length(lagTimes),length(fpasses));
% Cell for legend
legCell = {'Unfiltered', 'Linear fit subtracted'};
% Some constants
Kb = 1.38064852e-23; % Boltzmann constant
T = 273 + 20; % Assume trap at room temperature - maybe I should model this?
% radiusBead = 
% gamma = 6 .* pi .* radiusBead .* eta;
offset = length(legCell);

for fpIdx = 1:length(fpasses)
    fpass = fpasses(fpIdx);
    legCell{fpIdx + offset} = [num2str(fpass) 'Hz high-pass filtered'];
    %% Perform HP filtering
    % Crop time before HP
    % xCentresHP = highpass(xCentres(cropT(1):cropT(2)),fpass,fs) .* mPerPx;
    % yCentresHP = highpass(yCentres(cropT(1):cropT(2)),fpass,fs) .* mPerPx;
    
    % Crop time after HP
    xCentresHP = highpass( xCentres, fpass, fs) .* mPerPx;
    xCentresHP = xCentresHP(cropT(1):cropT(2));
    yCentresHP = highpass( yCentres, fpass, fs) .* mPerPx;
    yCentresHP = yCentresHP(cropT(1):cropT(2));
    
    xCentresHPcrop = xCentresHP(cropTHP(1):cropTHP(2));
    yCentresHPcrop = yCentresHP(cropTHP(1):cropTHP(2));
    
    %% Allan variancing
    
    [YHP(:,fpIdx), ~] = allanvar(xCentresHPcrop, lagTimes, fs);
end

[Y, Tau] = allanvar(xCentresM(cropTHP(1):cropTHP(2)), lagTimes, fs);
[Yuf, ~] = allanvar(xCentres(cropT(1):cropT(2)) .* mPerPx, lagTimes, fs);
% stdErr = sqrt((sqrt(2).* Kb .* T .*


fh = figure(20 + fileIdx);
fh.Name = fName;
clf

loglog(Tau, 1e9 .* sqrt([Yuf, Y, YHP]),'LineWidth',2)

title('Allan variance comparing high-pass to polynomial fit filtering')
legend(legCell, 'Location','best')
xlabel('Lag time, \tau (s)')
ylabel('Allan deviation (square root of Allan variance) (nm)')


