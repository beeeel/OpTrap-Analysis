function [fits, varargout] = unwrap_cell_v1(Imstack, centres, radii, varargin)
%UNWRAP_CELL_V1 - perform radial unwrapping of the cell using interpolation
%
%unwrap_cell_v1(Imstack, centres, radii, varargin) - if Imstack has
%N frames, centres and radii must be 2xN and 1xN and contain
%[X_centre, Y_centre] and radius in each column, respectively.
% 
% Draws n_theta different lines on each image space, from the centre given,
% to the maximum radius given, equally spaced around a full turn. Uses
% interp_method to perform interpolation, creating a polar-coordinates
% image matrix. Fits the polar equation for an ellipse, and uses
% post-processing to solve the equivalent equations fitting problem (a pi/2
% rotation and an inversion [b>a] produces the same curve, so half the fits
% are wrong by pi/2). Fitting is performed on data repeated sequentially to
% reduce fitting error (e.g.: if data is [1, 2, 3], fitting is performed on
% [1, 2, 3, 1, 2, 3, 1, 2, 3] when n_reps = 3)
%
% Returns a 3xN matrix of corrected fitted values [a; b; phi]. a and b are
% semi-major and semi-minor axis lengths, respectively. Optional output
% arguments are the unwrapped dataset and Ia, the filtered column maxima of
% the unwrapped dataset
% 
% Possible input arguments:
%   sc_up       - Scale factor for maximum sampling radius relative to maximum input radius
%   sc_down     - Scale factor for minimum sampling radius relative to maximum input radius
%   n_theta     - Number of angles to sample at (above 2*pi*r is oversampling)
%   n_reps      - Number of repeats for fitting data
%   tol         - Tolerance relative to median for discarding datapoints from fitting
%   inter_method- Interpolation method for unwrapping
%   centering   - Fit the off-centre equation instead of the centred one (experimental)
%   ifNaN       - What do if find_cell fails ('mean' or 'last')
%
%tol: For large deformation, a larg tol value is needed - for D=0.1, tol>=0.25
%
%ifNaN: when find_cell fails, unwrap_cell can either use the values for
% centre and radius from the previous frame (if frames from the start fail,
% they will not be unwrapped), or can use the average that find_cell did
% find.
%
%Note: This will throw "Error using reshape" if you try and run it with
% only one frame, or a single frame's worth of centres/radii.

fields = {'sc_up', 'n_theta', 'n_reps', 'tol', 'inter_method', 'sc_down',...
    'centering', 'ifNaN'};
defaults = {1.2, 360, 10, 0.15, 'linear', 0.5,...
    0, 'mean'};

par = cell2struct(defaults, fields,2);

% Parse inputs and create struct with parameters given
def_argin = 3;
if nargin > def_argin
    if mod(nargin,2) == 0 % Imstack counts towards nargin
        error('Please supply arguments in name-value pairs');
    end
    st = dbstack;
    % For each field, display it depending on its size and contents
    fprintf('Input arguments for %s:\n',st(1).name)
    for field = 1:(nargin - def_argin)/2
        if size(varargin{2 * field}, 2) < 4 && ...
                size(varargin{2 * field},1) == 1
            if isnumeric(varargin{2 * field})
                disp([varargin{2 * field - 1}, ' = ', ...
                    num2str(varargin{2 * field})]);
            elseif iscell(varargin{2 * field})
                fprintf([varargin{2 * field - 1}, ': '])
                for idx = 1:length(varargin{2 * field})/2
                    if isnumeric(varargin{2*field}{2*idx})
                        disp([varargin{2 * field}{2 * idx - 1}, ' = ', ...
                            num2str(varargin{2*field}{2*idx})]);
                    else
                        disp([varargin{2 * field}{2 * idx - 1}, ' = ', ...
                            varargin{2 * field}{2 * idx}])
                    end
                end
            end
        else
            disp([varargin{2 * field - 1}, ' = size[', ...
                num2str([size(varargin{2 * field},1), size(varargin{2 * field},2)])...
                , ']']);
        end
        % And put it into the par struct
        par.(varargin{2*field - 1}) = varargin{2*field};
    end
end

% Preallocate
r = (1:round(par.sc_up * max(radii)))';
fits = zeros(3 + 2 * par.centering,length(Imstack{1})); % get size correct if fitting for off-centre
theta = linspace(single(1),single(360),par.n_theta);

% Fix any NaNs from find_cell failing
if strcmp(par.ifNaN, 'last')
    % Take the last good value
    for frame = 1:length(radii)
        if ~isfinite(radii(frame)) && frame > 1
            radii(frame) = radii(frame - 1);
            centres(:, frame) = centres(:, frame - 1);
        end
    end
elseif strcmp(par.ifNaN, 'mean')
    % Take the mean value of those found 
    centres(:,~isfinite(radii)) = repmat(mean(centres(:,isfinite(radii)),2),1,sum(~isfinite(radii)));
    radii(~isfinite(radii)) = mean(radii(isfinite(radii)));
else
    error('ifNaN must have value ''mean'' or ''last''')
end

% From parametric equation of ellipse (x = a cos(t), y = b sin(t)), take to
% polar, r = (a.cos(f * theta + phi))^2 +  (b.sin(f * theta + phi))^2,
% where f is the frequency - twice per 2pi, theta is the angle (x below),
% and phi is the phase offset (orientation of ellipse). For an off-centre
% elipse, x = a cos(t) + dx, y = a sin(t) + dy.
if ~par.centering
    fiteqn = @(a, b, phi, x) sqrt((a * cos(0.0175*x + phi)).^2 + (b * sin(0.0175*x + phi)).^2 );
    lb = [0 0 0];
    ub = [inf inf pi];
    startval = [repmat(radii',1,2) repmat(1.5,size(radii,2),1)]; % Use radius from find_cell as a start point
elseif par.centering
    fiteqn = @(a, b, phi, dx, dy, x) sqrt((a * cos(0.0175*x + phi) + dx).^2 + (b * sin(0.0175*x + phi) + dy).^2 );
    lb = [0 0 0 -size(Imstack{1}{1,1})/2];
    ub = [inf inf pi size(Imstack{1}{1,1})/2];
    startval = [repmat(radii',1,2) repmat([1.5 0 0],size(radii,2),1)];
else
    error('Centering argument must have value 0 or 1')
end

tic
% For speed, everything is in one call to interp3, returning single which
% is then cast to uint16. The first argument is the whole imagestack, the
% remaining are x, y, z query points, all arguments are cast to single for
% memory efficiency
unwrapped = uint16(interp3(...
    single(cat(3,Imstack{1}{:,1})),...                                                              % V
    single(reshape(centres(1,:),1,1,length(centres)) + r.*cos(theta*pi/180)),...                    % Xq
    single(reshape(centres(2,:),1,1,length(centres)) + r.*sin(theta*pi/180)),...                    % Yq
    single(repmat(reshape(1:length(Imstack{1}),1,1,length(Imstack{1})),length(r),par.n_theta)),...  % Zq
    par.inter_method));                                                                             % METHOD
fprintf('Finished unwrapping in %gs\n',toc)

% Find the maxes and apply filtering based on deviation from median (will
% need large tolerance for large deformation) - maybe 0.25 for D = 0.1 
% To make this more robust, only check for maxes away from the centre (r=0)
% - do this by indexing a range scaled according to the "sc_down"
[~, Ia] = max(unwrapped(floor(par.sc_up*max(radii)*par.sc_down):end,:,:));
Ia = Ia + floor(par.sc_up*max(radii)*par.sc_down);
idxa = Ia > (1-par.tol)*median(Ia,2) & Ia < (1+par.tol)*median(Ia,2);
% For each frame, take angle corresponding to fit points and perform the
% fit
disp('Starting fitting')
for frame = 1:length(Imstack{1})
    % Take corresponding angles, repeat them and unwrap
    th_fit = repmat(theta(idxa(:,:,frame)),1,par.n_reps) + ...
        reshape(360*(0:par.n_reps-1).*ones(length(theta(idxa(:,:,frame))),1),1,[]);
    fitobj = fit(double(th_fit'),repmat(Ia(:,idxa(:,:,frame),frame)',par.n_reps,1),fiteqn,...
        'Lower', lb, 'Upper', ub, 'Start', startval(frame,:));
    fits(1:3,frame) = [fitobj.a, fitobj.b, fitobj.phi];
    if par.centering; fits(4:5,frame) = [fitobj.dx, fitobj.dy]; end
    prog = ceil(100 * frame / length(Imstack{1}));
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
end
% Fix the equivalent fitting equations problem - if b>a
fits(3,:) = fits(3,:) + (fits(2,:) > fits(1,:)) * pi/2;
fits(3,:) = mod(fits(3,:)+pi/2, pi) - pi/2;
fits(1:2,:) = [max(fits(1:2,:)); min(fits(1:2,:))].^2;
fprintf('%s\n',repmat(' ',1,104))
fprintf('Fitted data in %gs\n',toc)
switch nargout
    case 3; varargout = {unwrapped, Ia};
    case 4; varargout = {unwrapped, Ia, fiteqn};
end
end