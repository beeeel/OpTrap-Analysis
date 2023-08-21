function data = bead_preProcessCentres(data)
%% General preprocessing on bead data - polynomial fit removal, angular correction, unit conversion, time regularisation
% Takes options from data.opts. Options are
%   cropT   - crop time domain from index 1 to index 2
%   forceRun- force analysis code to run even if output field exists
%   pOrder  - polynomial fit order
%   angleCorrection - what it says (bead and cell data)
%   timeRegularisation - replace time vector with regular vector over same range
%   downsampleR - downsample position-time data by factor R
%   UseField - which field in data.pro to use (automatically updates)

if isfield(data,'mPerPx')
    mPerPx = data.mPerPx;
else
    mPerPx = 0.065e-6;
    data.mPerPx = mPerPx;
    warning('Using default value (%s m/px) for pixel size calibration', mPerPx)
end

if ~isfield(data, 'nPoints')
    data.nPoints = numel(data.raw.timeVecMs);
end

% Check opts exists before checking its contents
if ~isfield(data,'opts')
    data.opts = struct();
end
% Check essential fields have been put in opts
fnames = {'cropT', 'forceRun', 'pOrder', 'angleCorrection', ...
    'timeRegularisation', 'downsampleR', 'UseField','bandstop'};
defaults = {[1 data.nPoints], 0, 0, false, ...
    true, 1, '', []};
for fi = 1:length(fnames)
    if ~isfield(data.opts, fnames{fi})  || isempty(data.opts.(fnames{fi}))
        data.opts.(fnames{fi}) = defaults{fi};
    end
end

% We also should have an fName
if ~isfield(data,'fName')
    data.fName = '';
end

% Do time regularisation?
if data.opts.timeRegularisation
    t = data.raw.timeVecMs;
    dt = median(diff(t)) * data.opts.downsampleR;
    [~, I] = min(diff(t));
    % Gotta check for split acquisitions! These are characterised by a time
    % vector that resets to 0 partway, so this should do an okay job
    if sum(t==0) > 1
        error('Split acquisition detected! First acquisition ends at %i, second has length %i\ncd %s && datCropper.sh %i\t is a possible command to fix?', I, length(t)-I, data.dirPath,-I)
        % A useful command for trimming acquisitions: find . -type f -exec dd if={} of=../{} ibs=8000 count=523 \;
        % Note: First move data files into a new folder then run this with
        % either count or skip set as appropriate
    end
    data.pro.timeVecMs = ( ( 1:round(data.nPoints / data.opts.downsampleR ) ) - 1 ) * dt;
end

% These are in units of pixels
xCentres = data.raw.xCentresPx;
yCentres = data.raw.yCentresPx;

% Do angular correction?
if data.opts.angleCorrection
    N_doAngleCorrection;
elseif isfield(data, 'metadata')
    % Get ROI position
    roi = str2double( strsplit( data.metadata.FrameKey_0_0_0.ROI, '-' ) );
    data.opts.roi = roi;
end

% Store mean and std of co-ordinates AFTER angle correction, BEFORE demean.
data.pro.meanstd = [mean([xCentres; yCentres], 2) std([xCentres; yCentres], 0, 2)].* mPerPx;

if data.opts.pOrder > 0 && data.opts.angleCorrection
    warning('Be careful, drift removal was done ~before~ after conversion to angular co-ordinates,')
    warning('I have not thought carefully about the implications of this. Continue at own risk')
end
        
% Conditional drift removal only demeans when pOrder = 0
dims = [1, 3, 2];

Ts = data.opts.cropT(1):data.opts.cropT(2);

[~, xCentres, ~] = func_thermal_rm(Ts, ...
    permute(xCentres(:,Ts), dims), data.opts.pOrder, 1, length(Ts));
[~, yCentres, ~] = func_thermal_rm(Ts, ...
    permute(yCentres(:,Ts), dims), data.opts.pOrder, 1, length(Ts));

warning('I''ve changed how cropT is handled with data.pro... You''re gonna get errors in other bits of code, sorry!')

if data.opts.downsampleR > 1
    N_downsample
end

data.pro.xCentresM = ipermute(xCentres.* mPerPx, dims);
data.pro.yCentresM = ipermute(yCentres.* mPerPx, dims);
data.opts.UseField = 'CentresM';


function N_doAngleCorrection

    if isfield(data, 'ImstackFullFoV')
        % Get ROI position
        roi = str2double( strsplit( data.metadata.FrameKey_0_0_0.ROI, '-' ) );
        data.opts.roi = roi;
        
        % Get cell centre position
        cFile = [data.dirPath '/cell_centre.txt'];
        if exist(cFile, 'file')
            % Load previously measured cell centre
            fprintf('Loading cell centre from file: %s\n', cFile)
            cCentre = load(cFile, '-ascii');
        else
            % Measure centre function
            cCentre = measure_cell_centre(data.ImstackFullFoV{1}{1,1}, data.dirPath, roi);
        end
        data.opts.cCentre = cCentre;
        % This is written in my lab book 22/9/2021
        rhoX = xCentres + roi(1) - cCentre(1);
        rhoY = yCentres + roi(2) - cCentre(2);
        % This is actually radial co-ordinate
        xCentres = sqrt( rhoX.^2 + rhoY.^2 );
        % This is r times tangential co-ordinate
        yCentres = xCentres .* unwrap(2*atan(rhoY ./ rhoX))/2; % Gotta double and then half to make unwrap work.
    else
        data.opts.angleCorrection = false;
        warning('No image files loaded, skipping angular correction')
    end
end

    function N_downsample
        % Downsample position and time data, storing the outputs. Untested.
        [nR, ~, nT] = size(xCentres);
        r = data.opts.downsampleR;
        
        if all(data.opts.cropT == [1 nT])
            cropT = round([1 nT./r]);
        else
            cropT = round(data.opts.cropT./r);
            if cropT(1) == 0
                cropT(1) = 1;
            end
        end
        data.opts.cropT = cropT;
            
        validateattributes(r, {'numeric'}, {'integer','positive'})
        
        tmp = [xCentres, yCentres];
        Centres = zeros(nR, 2, floor(nT./r));
        for idx = 1:floor(nT/r)
            % Average both rows of the track, using a number of elements
            % equal to R (downsampling factor)
            Centres(:,:,idx) = mean( tmp( :, :, (idx - 1) * r + 1 : idx * r ),3);
        end
        xCentres = Centres(:,1,:);
        yCentres = Centres(:,2,:);
        
        if ~data.opts.timeRegularisation
            tmp = data.pro.timeVecMs;
            t1 = zeros(1, floor(nT./r));
            for idx = 1:floor(nT/r)
                % Average time vector using a number of elements equal to R
                % (downsampling factor)
                t1(idx) = mean( tmp( (idx - 1) * r + 1 : idx * r ));
            end
            data.pro.timeVecMs = t1(cropT(1):cropT(2));
            data.raw.timeVecMs = t1;
        end
        
        data.opts.cropT = [1 floor(nT./r)];
    end
end
