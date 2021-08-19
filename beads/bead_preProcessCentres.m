function data = bead_preProcessCentres(data)

if isfield(data,'mPerPx')
    mPerPx = data.mPerPx;
else
    warning('Using default value for pixel size calibration')
    mPerPx = 0.07e-6;
end

% Do angular correction?
if ~isfield(data.opts, 'angleCorrection')
    data.opts.angleCorrection = false;
end

xCentresM = data.raw.xCentresPx .* mPerPx;
yCentresM = data.raw.yCentresPx .* mPerPx;

% Conditional drift removal only demeans when pOrder = 0
dims = [1, 3, 2];

[~, xCentresM, ~] = func_thermal_rm(1:length(xCentresM), ...
    permute(xCentresM, dims), data.opts.pOrder, 1, length(xCentresM));
[~, yCentresM, ~] = func_thermal_rm(1:length(yCentresM), ...
    permute(yCentresM, dims), data.opts.pOrder, 1, length(yCentresM));

if data.opts.angleCorrection
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
        % Apply correction for ROI position
        roi = str2double( strsplit( data.metadata.FrameKey_0_0_0.ROI, '-' ) );
        cCentre = cCentre - roi(1:2);
        % This is actually radial co-ordinate
        xCentresM = sqrt( ( xCentresM + cCentre(1) ).^2 ...
            + ( yCentresM + cCentre(2) ).^2 );
        % This is r times tangential co-ordinate
        yCentresM = xCentresM .* atan( ( yCentresM + cCentre(2) ) ...
            ./ ( xCentresM + cCentre(1) ) );
        if data.opts.pOrder > 0
            warning('Be careful, drift removal was done before conversion to angular co-ordinates,')
            warning('I have not thought carefully about the implications of this. Continue at own risk')
        end
    else
        warning('No image files loaded, skipping angular correction')
    end
end

data.pro.xCentresM = ipermute(xCentresM, dims);
data.pro.yCentresM = ipermute(yCentresM, dims);
data.opts.UseField = 'CentresM';