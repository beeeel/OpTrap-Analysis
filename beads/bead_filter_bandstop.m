function data = bead_filter_bandstop(data, doPlots)
%% data = bead_filter_bandstop(data)
% Performs bandstop filtering of data. Use data.opts.bandstop to define 1
% band per row ([FpassLower FpassUpper]).

%% Setup

if ~isfield(data.opts, 'bandstop')
    error('Need field bandstop in  data.opts')
elseif size(data.opts.bandstop, 2) ~= 2
    error('Expected width of data.opts.bandstop to be 2, instead found %i\n', size(data.opts.bandstop, 2))
end

autoBand = isfield(data.opts, 'autoband');

if ~exist('doPlots', 'var')
    doPlots = true;
end

if ~isfield(data.opts, 'centresRow')
    centresRow = 1;
else
    centresRow = data.opts.centresRow;
end

mPerPx = 65e-9;
if isfield(data, 'mPerPx')
    mPerPx = data.mPerPx;
end

cropT = [1 data.nPoints];
if isfield(data.opts, 'cropT')
    if ~isempty(data.opts.cropT)
        cropT = data.opts.cropT;
    end
end

% Take 1 row of centres and sandwich it with mirrored versions to reduce
% artefacts is commented out
if ~isfield(data.opts, 'UseField')
    tmp = mPerPx * ...[data.raw.xCentresPx(centresRow, cropT(2):-1:cropT(1)) ...
        data.raw.xCentresPx(centresRow, cropT(1):cropT(2)) ...
        ;% data.raw.xCentresPx(centresRow, cropT(2):-1:cropT(1))];
    centres(:,1) = tmp;
    
    tmp = mPerPx * ...[data.raw.yCentresPx(centresRow, cropT(2):-1:cropT(1)) ...
        data.raw.yCentresPx(centresRow, cropT(1):cropT(2)) ...
        ;%data.raw.yCentresPx(centresRow, cropT(2):-1:cropT(1))];
    centres(:,2) = tmp;
else
    fN = data.opts.UseField;
    tmp = mPerPx * ...[data.pro.(['x' fN ])(centresRow, cropT(2):-1:cropT(1)) ...
        data.pro.(['x' fN ])(centresRow, cropT(1):cropT(2)) ...
        ;%data.pro.(['x' fN ])(centresRow, cropT(2):-1:cropT(1))];
    centres(:,1) = tmp;
    
    tmp = mPerPx * ...[data.pro.(['y' fN ])(centresRow, cropT(2):-1:cropT(1)) ...
        data.pro.(['y' fN ])(centresRow, cropT(1):cropT(2)) ...
        ;%data.pro.(['y' fN ])(centresRow, cropT(2):-1:cropT(1))];
    centres(:,2) = tmp;
end

fPasses = data.opts.bandstop;
fS = 1/diff(data.raw.timeVecMs(1:2)*1e-3);

for fPassIdx = 1:size(fPasses,1)
    if ~autoBand
        [centres, ~] = bandstop(centres, fPasses(fPassIdx, :), fS);
    else
        error('Sorry autoBand has not been programmed yet')
    end
end

% centres = centres(cropT(2)+1:2*cropT(2),:);

data.pro.xCentresBS = centres(:,1)';
data.pro.yCentresBS = centres(:,2)';

data.opts.UseField = 'CentresBS';

if doPlots
    t = data.raw.timeVecMs(cropT(1):cropT(2))*1e-3;
    xy = 'XY';
    
    figure
    
    for plt = 1:2
        subplot(2, 2, plt);
        plot(t, centres(:,plt)*1e6);
        title(sprintf('%s centres time trace', xy(plt)))
        xlabel('Time (s)')
        ylabel('Centre position (μm)')
    end
    
    for plt = 1:2
        ax = subplot(2, 2, 2 + plt);
        fft_scaled(t, 1e6 * centres(:, plt)', true, ax);
        ax.XLim(1) = 10;
        ylabel('Amplitude (nm)')
        drawnow
    end
    
end
end