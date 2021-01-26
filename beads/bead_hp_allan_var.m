function data = bead_hp_allan_var(data, field, fpasses, cropTHPval, doPlot)
%% data = bead_hp_allan_var(data, field, fpasses, cropTHPval, doPlot)
%% High-pass filter position and then calculate allan variance. 
% Do this for each frequency in fpasses and (optional) plot all the
% variances together, along with unfiltered and linear fitted.
%
% Uses position stored in data.raw.(field)


centreVec = data.raw.(field);
timeVec = data.raw.timeVecMs;
cropT = data.opts.cropT;
fName = data.fName;

% Outputs

% Need to change input validator
%% Setup
fs = length(centreVec)./ (1e-3 * (max(timeVec) - min(timeVec))); % Sampling frequency in Hz
%fpasses = [0.5 1 5 10 50 100]; % Pass frequency
%cropTHPval = 1; % Number of frames to crop to remove ringing caused by HP

%% Input validation
N_input_validator();

%% Preallocate
% Crop indices for after high-passing (remove ringing)
cropTHP = [cropTHPval+1, diff(cropT) + 1 - cropTHPval];
% Number of points in final HP vector
nPoints = diff(cropT) - 2 * cropTHPval + 2;
% Lag times
lagTimes = calcLagTimes(nPoints);
YOut = zeros(length(lagTimes),length(fpasses));
% Cell for legend
legCell = {'Unfiltered', 'Linear fit subtracted'};
legOffset = length(legCell);
if isfield(data,'mPerPx')
    mPerPx = data.mPerPx;
else
    warning('Using default value for pixel size calibration')
    mPerPx = 0.07e-6;
end

for fpIdx = 1:length(fpasses)
    fpass = fpasses(fpIdx);
    legCell{fpIdx + legOffset} = [num2str(fpass) 'Hz high-pass filtered'];
    %% Perform HP filtering
    % Crop time before HP
    % xCentresHP = highpass(xCentres(cropT(1):cropT(2)),fpass,fs) .* mPerPx;
    % yCentresHP = highpass(yCentres(cropT(1):cropT(2)),fpass,fs) .* mPerPx;
    
    % Crop time (twice) after HP
    CentresHP = highpass( centreVec, fpass, fs) .* mPerPx;
    CentresHP = CentresHP(cropT(1):cropT(2));
    CentresHPcrop = CentresHP(cropTHP(1):cropTHP(2));
    
    %% Allan variancing
    
    [YOut(:,fpIdx), ~] = allanvar(CentresHPcrop, lagTimes, fs);
end
xCentresM = centreVec(cropT(1):cropT(2)) .* mPerPx;
dims = [1, 3, 2];
pOrder = 1;
[~, xCentresM, ~] = func_thermal_rm(1:length(xCentresM), ...
    permute(xCentresM, dims), pOrder, 1, length(xCentresM));

xCentresM = ipermute(xCentresM, dims);

[Y, Tau] = allanvar(xCentresM(cropTHP(1):cropTHP(2)), lagTimes, fs);
[Yuf, ~] = allanvar(centreVec(cropT(1):cropT(2)) .* mPerPx, lagTimes, fs);

data.pro.([field(1) 'AlanVar']) = YOut;
data.pro.([field(1) 'Tau']) = Tau;
data.pro.([field(1) 'CentresHP']) = CentresHPcrop;


if doPlot
    fh = figure;
    fh.Name = fName;
    clf
    
    loglog(Tau, 1e18 .* [Yuf, Y, YOut],'LineWidth',2)
    
    title('Allan variance comparing high-pass to polynomial fit filtering')
    legend(legCell, 'Location','best')
    xlabel('Lag time, \tau (s)')
    ylabel('Allan variance (nm^2)')
    drawnow
end


function N_input_validator()
% Inputs are: (centreVec, timeVec, fpasses, cropT, cropTHPval, fName, doPlot)
validateattributes(centreVec, {'numeric'},{'vector','real','finite'},...
    'func_bead_hp_allan_var','centreVec')
validateattributes(timeVec, {'numeric'},...
    {'nonnegative','increasing','vector','real','finite','size',size(centreVec)},...
    'func_bead_hp_allan_var','timeVec')
validateattributes(fpasses, {'numeric'},{'vector','real','finite', 'positive','<',fs},...
    'func_bead_hp_allan_var','fpasses')
validateattributes(cropT, {'numeric'},{'numel',2,'integer', 'positive','<=',length(centreVec)},...
    'func_bead_hp_allan_var','cropT')
validateattributes(cropTHPval, {'numeric'},{'integer','scalar', 'positive','<',(diff(cropT)+1)/2},...
    'func_bead_hp_allan_var','cropTHPval')
validateattributes(fName,{'string','char'},{},'func_bead_hp_allan_var','fName')
validateattributes(doPlot,{'logical'},{'scalar'},'func_bead_hp_allan_var','doPlot')
end

end