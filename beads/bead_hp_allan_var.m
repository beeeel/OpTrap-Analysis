function data = bead_hp_allan_var(data, field, fpasses, cropTHPval, doPlot, varargin)
%% data = bead_hp_allan_var(data, field, fpasses, cropTHPval, doPlot, [centresRow])
%% High-pass filter position and then calculate allan variance. 
% Do this for each frequency in fpasses and (optional) plot all the
% variances together, along with unfiltered and linear fitted.
%
% Uses position stored in data.raw.(field)
% (Lazily) uses first row of centres

error('Check how cropTHP is handled because you''re phasing it out')

% Choose centres row
if isfield(data.opts, 'centresRow')
    cRow = data.opts.centresRow;
else
    cRow = 1;
end
% User input overrides for legacy porpoises
if nargin >= 6
    cRow = varargin{1};
end

data.opts.([field(1) 'HPSuffix']) = data.raw.suffixes(cRow);

if isfield(data.raw, field)
    centreVec = data.raw.(field);
    % If we're using raw, we need to apply calibration
    if isfield(data,'mPerPx')
        mPerPx = data.mPerPx;
    else
        warning('Using default value for pixel size calibration')
        mPerPx = 0.065e-6;
    end
    
    timeVec = data.raw.timeVecMs;
    
elseif isfield(data.pro, field)
    centreVec = data.pro.(field);
    % If we're using pro, we mustn't reapply calibration
    mPerPx = 1;
    
    timeVec = data.pro.timeVecMs;

else
    error('Could not find field %s in data.raw or data.pro',field)
end

centreVec = [centreVec(cRow,end:-1:1) centreVec(cRow,:) centreVec(cRow,end:-1:1)];


if length(data.opts.cropT) == 2
    cropT = data.opts.cropT;      
else
    cropT = [1 length(timeVec)];
end

cropT = cropT + length(timeVec);
fName = data.fName;

% Outputs

% Need to change input validator
%% Setup
fs = length(timeVec)./ (1e-3 * (max(timeVec) - min(timeVec))); % Sampling frequency in Hz

%% Input validation
N_input_validator();

%% Preallocate
% Crop indices for after high-passing (remove ringing)
cropTHP = [cropTHPval+1, diff(cropT) + 1 - cropTHPval];

% % Number of points in final HP vector
% %nPoints = diff(cropT) - 2 * cropTHPval + 2;
% % Lag times
% %lagTimes = calcLagTimes(nPoints);
% %YOut = zeros(length(lagTimes), size(centreVec,1), length(fpasses));

% Cell for legend
legCell = {'Unfiltered', 'Linear fit subtracted'};
legOffset = length(legCell);

for fpIdx = 1:length(fpasses)
    fpass = fpasses(fpIdx);
    legCell{fpIdx + legOffset} = [num2str(fpass) 'Hz high-pass filtered'];
    %% Perform HP filtering
    % Crop time before HP
    % xCentresHP = highpass(xCentres(cropT(1):cropT(2)),fpass,fs) .* mPerPx;
    % yCentresHP = highpass(yCentres(cropT(1):cropT(2)),fpass,fs) .* mPerPx;
    
    % Crop time (twice) after HP
    CentresHP = highpass( centreVec', fpass, fs)' .* mPerPx;
    CentresHP = CentresHP(:,cropT(1):cropT(2));
    CentresHPcrop = CentresHP(:,cropTHP(1):cropTHP(2));
    
    %% Allan variancing
    
%    [YOut(:,:,fpIdx), Tau] = allanvar(CentresHPcrop', lagTimes, fs);
end


%data.pro.([field(1) 'AlanVar']) = YOut;
%data.pro.([field(1) 'Tau']) = Tau;
data.pro.([field(1) 'CentresHP']) = CentresHPcrop;
data.opts.UseField = 'CentresHP';

if doPlot
    warning('skipping plot to avoid calling allanvar (because it''s not installed on string')
end
if false
    xCentresM = centreVec(cropT(1):cropT(2)) .* mPerPx;
    dims = [1, 3, 2];
    pOrder = 1;
    [~, xCentresM, ~] = func_thermal_rm(1:length(xCentresM), ...
        permute(xCentresM, dims), pOrder, 1, length(xCentresM));
    
    xCentresM = ipermute(xCentresM, dims);
    
    [Y, ~] = allanvar(xCentresM(cropTHP(1):cropTHP(2)), lagTimes, fs);
    [Yuf, ~] = allanvar(centreVec(cropT(1):cropT(2)) .* mPerPx, lagTimes, fs);

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
validateattributes(centreVec, {'numeric'},{'real','finite'},...
    'func_bead_hp_allan_var','centreVec')
validateattributes(timeVec, {'numeric'},...
    {'nonnegative','increasing','vector','real','finite'},...,'size',size(centreVec)},...
    'func_bead_hp_allan_var','timeVec')
validateattributes(fpasses, {'numeric'},{'vector','real','finite', 'positive','<',fs},...
    'func_bead_hp_allan_var','fpasses')
validateattributes(cropT, {'numeric'},{'numel',2,'integer', 'positive','<=',length(centreVec)},...
    'func_bead_hp_allan_var','cropT')
validateattributes(cropTHPval, {'numeric'},{'integer','scalar', 'nonnegative','<',(diff(cropT)+1)/2},...
    'func_bead_hp_allan_var','cropTHPval')
validateattributes(fName,{'string','char'},{},'func_bead_hp_allan_var','fName')
validateattributes(doPlot,{'logical'},{'scalar'},'func_bead_hp_allan_var','doPlot')
end

function lT = calcLagTimes(nP, varargin)
%% Calculate lag times to use with allanvar
% lT = calcLagTimes(nPoints, [n])

if nargin == 1
    n = 128;
elseif nargin == 2
    n = varargin{1};
else
    error('Too many nargins!')
end

% This next line is a shitshow. I want logarithmically spaced points
% between 1 and the half length of my centres vector. Logspace wants the
% decades (i.e.: 5e5 is in decade 5), so use log10 to get that. Allanvar
% wants a list of unique integers, so round and remove duplicates.
lT = unique(round(logspace(0,floor(log10(nP/2)),n)));
end

end
