function data = bead_preProcessCentres(data)
%% General preprocessing on bead data - polynomial fit removal, angular correction, unit conversion, time regularisation
if isfield(data,'mPerPx')
    mPerPx = data.mPerPx;
else
    warning('Using default value for pixel size calibration')
    mPerPx = 0.07e-6;
end

% Do time regularisation?
if ~isfield(data.opts, 'timeRegularisation')
    data.opts.timeRegularisation = false;
elseif data.opts.timeRegularisation
    t = data.raw.timeVecMs;
    dt = median(diff(t));
    data.raw.timeVecMs = ( (1:data.nPoints) - 1 ) * dt;
end

% These are in units of pixels
xCentres = data.raw.xCentresPx;
yCentres = data.raw.yCentresPx;

% Do angular correction?
if ~isfield(data.opts, 'angleCorrection')
    data.opts.angleCorrection = false;
elseif data.opts.angleCorrection
    N_doAngleCorrection;
end

if data.opts.pOrder > 0 && data.opts.angleCorrection
    warning('Be careful, drift removal was done ~before~ after conversion to angular co-ordinates,')
    warning('I have not thought carefully about the implications of this. Continue at own risk')
end
        
% Conditional drift removal only demeans when pOrder = 0
dims = [1, 3, 2];

[~, xCentres, ~] = func_thermal_rm(1:length(xCentres), ...
    permute(xCentres, dims), data.opts.pOrder, 1, length(xCentres));
[~, yCentres, ~] = func_thermal_rm(1:length(yCentres), ...
    permute(yCentres, dims), data.opts.pOrder, 1, length(yCentres));

if ~isfield(data.opts.downsampleR)
    data.opts.downsampleR = 1;
elseif data.opts.downsampleR > 1
    N_downsample
end

data.pro.xCentresM = ipermute(xCentres.* mPerPx, dims);
data.pro.yCentresM = ipermute(yCentres.* mPerPx, dims);
data.opts.UseField = 'CentresM';


function N_doAngleCorrection

    if isfield(data, 'ImstackFullFoV')
        cFile = [data.dirPath '/cell_centre.txt'];
        if exist(cFile, 'file')
            % Load previously measured cell centre
            fprintf('Loading cell centre from file: %s\n', cFile)
            cCentre = load(cFile, '-ascii');
        else
            % Measure centre function
            cCentre = measure_cell_centre(data.ImstackFullFoV{1}{1,1}, data.dirPath);
        end
        % Get ROI position
        roi = str2double( strsplit( data.metadata.FrameKey_0_0_0.ROI, '-' ) );
        % This is written in my lab book 22/9/2021
        rhoX = xCentres + roi(1) - cCentre(1);
        rhoY = yCentres + roi(2) - cCentre(2);
        % This is actually radial co-ordinate
        xCentres = sqrt( rhoX.^2 + rhoY.^2 );
        % This is r times tangential co-ordinate
        yCentres = xCentres .* atan( rhoY ./ rhoX );
    else
        data.opts.angleCorrection = false;
        warning('No image files loaded, skipping angular correction')
    end
end

    function N_downsample
        % Downsample position and time data, storing the outputs. Untested.
        
        [nR, ~, nT] = size(xCentres);
        r = data.opts.downsampleR;
        
        tmp = [xCentres, yCentres];
        Centres = zeros(nR, 2, floor(nT./r));
        for idx = 1:floor(nT/r)
            % Average both rows of the track, using a number of elements
            % equal to R (downsampling factor)
            Centres(:,:,idx) = mean( tmp( :, :, (idx - 1) * r + 1 : idx * r ));
        end
        xCentres = Centres(:,1,:);
        yCentres = Centres(:,2,:);
        
        tmp = data.raw.timeVecMs;
        t1 = zeros(1, floor(nT./r));
        for idx = 1:floor(nT/r)
            % Average time vector using a number of elements equal to R
            % (downsampling factor)
            t1(idx) = mean( tmp( (idx - 1) * r + 1 : idx * r ));
        end
        data.raw.timeVecMs = t1;
    end
end