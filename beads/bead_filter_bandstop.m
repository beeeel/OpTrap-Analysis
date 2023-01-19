function data = bead_filter_bandstop(data, doPlots)
%% data = bead_filter_bandstop(data)
% Performs bandstop filtering of data. Use data.opts.bandstop to define 1
% band per row ([FpassLower FpassUpper]).

%% Setup

if ~isfield(data.opts, 'bandstop')
    error('Need field bandstop in  data.opts')
elseif isempty(data.opts.bandstop)
    warning('Skipping bandstop filtering')
    return
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

for rIdx = centresRow
    % Take 1 row of centres and sandwich it with mirrored versions to reduce
    % artefacts is commented out
    if ~isfield(data.opts, 'UseField')
        tmp = mPerPx * ...[data.raw.xCentresPx(centresRow, cropT(2):-1:cropT(1)) ...
            data.raw.xCentresPx(rIdx, cropT(1):cropT(2)) ...
            ;% data.raw.xCentresPx(centresRow, cropT(2):-1:cropT(1))];
        centres(:,1) = tmp;
        
        tmp = mPerPx * ...[data.raw.yCentresPx(centresRow, cropT(2):-1:cropT(1)) ...
            data.raw.yCentresPx(rIdx, cropT(1):cropT(2)) ...
            ;%data.raw.yCentresPx(centresRow, cropT(2):-1:cropT(1))];
        centres(:,2) = tmp;
        
        t = data.raw.timeVecMs*1e-3;
    else
        fN = data.opts.UseField;
        % Gotta assume that pixel calibration has already been applied - after
        % all we're taking the data from pro!
        tmp = data.pro.(['x' fN ])(rIdx, cropT(1):cropT(2));
        centres(:,1) = tmp;
        
        tmp = data.pro.(['y' fN ])(rIdx, cropT(1):cropT(2));
        centres(:,2) = tmp;
        
        t = data.pro.timeVecMs*1e-3;
    end
    
    fPasses = data.opts.bandstop;
    fS = 1/diff(t(1:2));
    if any(isnan(centres),'all')
        warning('Skipping bandstop filtering for field %s row %i\n',fN,rIdx)
    else
        for fPassIdx = 1:size(fPasses,1)
            if ~autoBand
                [centres, ~] = bandstop(centres, fPasses(fPassIdx, :), fS);
            else
                error('Sorry autoBand has not been programmed yet')
            end
        end
    end
    
    % centres = centres(cropT(2)+1:2*cropT(2),:);
    
    % This should put unfiltered data if there are NaNs
    data.pro.xCentresBS(rIdx,:) = centres(:,1)';
    data.pro.yCentresBS(rIdx,:) = centres(:,2)';
end

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