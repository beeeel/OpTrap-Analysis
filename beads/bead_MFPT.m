function data = bead_MFPT(data, nD, doPlots)

narginchk(1,3)

% If we're not given number of distances to calculate, use default
if ~exist('nD', 'var')
    nD = 100;
end

% If we're told to plot, don't
if ~exist('doPlots', 'var')
    doPlots = false;
end

cropT = data.opts.cropT;

% If there's a field we should use, use it
if isfield(data.opts, 'UseField')
    fn = data.opts.UseField;
elseif isfield(data.pro, 'xCentresM')
    fn = 'CentresM';
end

% If working with 1bead data, use 1 row, for 2bead data, take 2.
centresRow = 1;
if strcmp(data.raw.suffixes{1}, 'l') && strcmp(data.raw.suffixes{2}, 'r')
    centresRow = [1 2];
end

x = data.pro.(['x' fn])(centresRow,:);
y = data.pro.(['y' fn])(centresRow,:);
t = data.pro.timeVecMs*1e-3;

dMin = 1e-9;
dMax = range(x);

xFPi = nan(nD,size(x,2)-1);
% yFPT = zeros(size(x,2)-1,size(Ds,2));
% nFPT = 0;
% nt = size(x,2);

% Not 100% sure about this: 
% For each initial observation
parfor idx = 1:size(x,2) - 1
    % Create local variable for distances
    Ds = logspace(floor(log10(dMin)), ceil(log10(dMax)), nD);
    % Reset continuation flag
    cont = true;
    % Initialise first passage index array
    fpi = nan(nD,1);
    % Start with first subsequent observation
    jdx = idx + 1;
    % and first distance
    kdx = 1;
    D = Ds(kdx);
    % Keep going through subsequent observations
    while jdx < size(x,2) && cont
        % Calculating displacement
        dx = abs(x(idx) - x(jdx));
        % If the displacement is bigger than the passage distance
        while dx > D && cont
            % Record the first passage indices
            fpi(kdx) = (jdx - idx);
            % Increment the passage distance
            kdx = kdx + 1;
            cont = kdx <= nD;
            D = Ds(kdx);
            % This will end loop when D > dx or we've done all the passage
            % distances
        end
        % Then it's time to increment the subsequent observation index
        jdx = jdx + 1;
    end
    % Return complete first passage index array
    xFPi(:,idx) = fpi;
end

xMFPT = mean(xFPi.*dt,2,'omitnan');
