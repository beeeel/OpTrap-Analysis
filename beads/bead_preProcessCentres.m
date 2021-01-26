function data = bead_preProcessCentres(data)

if isfield(data,'mPerPx')
    mPerPx = data.mPerPx;
else
    warning('Using default value for pixel size calibration')
    mPerPx = 0.07e-6;
end

xCentresM = data.raw.xCentresPx(data.opts.cropT(1):data.opts.cropT(2)) .* mPerPx;
yCentresM = data.raw.yCentresPx(data.opts.cropT(1):data.opts.cropT(2)) .* mPerPx;

% Conditional drift removal only demeans when pOrder = 0
dims = [1, 3, 2];

[~, xCentresM, ~] = func_thermal_rm(1:length(xCentresM), ...
    permute(xCentresM, dims), data.opts.pOrder, 1, length(xCentresM));
[~, yCentresM, ~] = func_thermal_rm(1:length(yCentresM), ...
    permute(yCentresM, dims), data.opts.pOrder, 1, length(yCentresM));
data.pro.xCentresM = ipermute(xCentresM, dims);
data.pro.yCentresM = ipermute(yCentresM, dims);