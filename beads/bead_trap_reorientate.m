function data = bead_trap_reorientate(data, varargin)
%% data = bead_trap_reorientate(data, [doPlot])
% Perform co-ordinate transform to put "x" along weakest trap direction
% Use histogram of trap centres and regionprops to find trap anisotropy


if isfield(data.opts, 'UseField')
    fn = {'pro', data.opts.UseField};
    sf = 1;
elseif isfield(data.pro, 'CentresM')
    fn = {'pro', 'CentresM'};
    sf = 1;
else
    fn = {'raw', 'CentresPx'};
    if isfield(data, 'mPerPx')
        sf = data.mPerPx;
    else
        sf = 6.5e-8;
        warning('Using default pixel calibration 0.065 px/μm')
    end
end

if isfield(data.opts, 'centresRow')
    cR = data.opts.centresRow;
else
    cR = 1;
end

if nargin > 1
    doPlot = varargin{1};
else
    doPlot = true;
end

% Get the data
X = data.(fn{1}).(['x' fn{2}])(cR, :) .* sf;
Y = data.(fn{1}).(['y' fn{2}])(cR, :) .* sf;

% Do the histogramming
[N, xE, yE] = histcounts2(X, Y);

% Find a threshold - ignore least occupied 5% of bins
Nsort = sort(N(:));
Ncum = cumsum(Nsort)./sum(N(:));
idx = find(Ncum > 0.05,1);
binIm = N>Nsort(idx);

% Fix holes by dilate and erode
se = strel('disk', 2);
binIm = imdilate(binIm, se);
binIm = imerode (binIm, se);

% Get details from regionprops
stats = regionprops(binIm, N, {'Eccentricity', 'MajoraxisLength','MinoraxisLength','Orientation'});
data.pro.anisotropyStats = stats;

% Do the transformation
Xr = cosd(stats.Orientation) .* X - sind(stats.Orientation) .* Y;
Yr = +sind(stats.Orientation) .* X + cosd(stats.Orientation) .* Y;

data.pro.xCentresMr = Xr;
data.pro.yCentresMr = Yr;

if doPlot
    figure
    histogram2(Xr, Yr);
    xlabel('X (m)')
    ylabel('Y (m)')
end