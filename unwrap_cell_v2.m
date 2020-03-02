function [Fits, varargout] = unwrap_cell_v2(Imstack, Centres, Radius, varargin)
%UNWRAP_CELL_V2 - perform radial unwrapping of the cell using
%interpolation, handle off-centre by fitting twice.
% 
%unwrap_cell_v2(Imstack, centres, radii, varargin) - if Imstack has
%N frames, centres and radii must be 2xN and 1xN and contain
%[X_centre, Y_centre] and radius in each column, respectively.
%
% Differs from V1 in that it first fits for an off-centre circle, then for
% a centred ellipse.
% 
% Draws n_theta different lines on each image space, from the centre given,
% to the maximum radius given, equally spaced around a full turn. Uses
% interp_method to perform interpolation, creating a polar-coordinates
% image matrix. First fits for a circle, then fits the polar equation for
% an ellipse, and uses post-processing to solve the equivalent equations
% fitting problem (a pi/2 rotation and an inversion [b>a] produces the same
% curve, so half the fits are wrong by pi/2). Fitting is performed on data
% repeated sequentially to reduce fitting error (e.g.: if data is [1, 2,
% 3], fitting is performed on [1, 2, 3, 1, 2, 3, 1, 2, 3] when n_reps = 3)
%
% Returns a 3xN matrix of corrected fitted values [a; b; phi]. a and b are
% semi-major and semi-minor axis lengths, respectively. Additional optional
% output arguments below:
% 
% Varargout: 2 to 4 additional arguments, in this order: {Unwrapped, Ia,
% FitEqn, Offset}. Unwrapped is size [n_r, n_theta, n_frames] containing
% unwrapped cell images, Ia is size [1, n_theta, n_frames] containing
% indexes used for fitting, FitEqn is a function handle used for fitting
% the ellipse, Offset is size [3, n_frames], containing [r; dx; dy] where r
% is fitted circle radius, and dx, dy are fitted centre offsets
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
defaults = {1.2, 360, 5, 0.15, 'linear', 0.5,...
    0, 'mean'};

Par = cell2struct(defaults, fields,2);
def_argin = 3;

% Parse inputs and create struct with parameters given
if ~isempty(varargin)
    if mod(nargin,2) == 0 % Imstack counts towards nargin
        error('Please supply arguments in name-value pairs');
    end
    stack = dbstack;
    % For each field, display it depending on its size and contents
    fprintf('Input arguments for %s:\n',stack(1).name)
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
        Par.(varargin{2*field - 1}) = varargin{2*field};
    end
end

% Fix any NaNs from find_cell failing
switch Par.ifNaN
    case 'last' % Take the last good value
        for frame = 1:length(Radius)
            if ~isfinite(Radius(frame)) && frame > 1
                Radius(frame) = Radius(frame - 1);
                Centres(:, frame) = Centres(:, frame - 1);
            end
        end
    case 'mean'
        % Take the mean value of those found
        Centres(:,~isfinite(Radius)) = repmat(mean(Centres(:,isfinite(Radius)),2),1,sum(~isfinite(Radius)));
        Radius(~isfinite(Radius)) = mean(Radius(isfinite(Radius)));
    case 'centre'
        % Assume the cell is centred and fills the FoV
        Centres(:,~isfinite(Radius)) = repmat([ImH/2 ImW/2],1,sum(~isfinite(Radius)));
        Radius(~isfinite(Radius)) = (ImW < ImH) * ImW/2 + (ImW > ImH) * ImH/2;
    otherwise
        error('ifNaN must have value ''mean'' or ''last'' or ''centre''')
end

tic

% Parametric eqn of off-centre circle (x = cos(t) + dx, y = sin(t) + dy),
% convert to polar (r = x.^2 + y.^2)
CircleEqn = @(r, dx, dy, x) sqrt((r .* cos(0.0175*x) + dx).^2 + (r .* sin(0.0175*x) + dy).^2);
lb = [0, -size(Imstack{1}{1,1})/2];
ub = [Par.sc_up * max(Radius), size(Imstack{1}{1,1})/2];
StartVal = [Radius', repmat([0, 0],size(Radius,2), 1)];

% Perform unwrapping and fitting for off-centreness
[Offset, ~, ~] = UnwrapAndFit(Imstack, CircleEqn, Radius,  Centres, lb, ub, StartVal, Par);

% From parametric equation of ellipse (x = a cos(t), y = b sin(t)), take to
% polar, r = (a.cos(f * theta + phi))^2 +  (b.sin(f * theta + phi))^2,
% where f is the frequency - once per 2pi (doubled by squaring), theta is
% the angle (x below is in degrees), and phi is the phase offset
% (orientation of ellipse). For an off-centre elipse, x = a cos(t) + dx, y
% = a sin(t) + dy.
if ~Par.centering
    FitEqn = @(a, b, phi, x) sqrt((a * cos(0.0175*x + phi)).^2 + (b * sin(0.0175*x + phi)).^2 );
    lb = [0 0 0];
    ub = [inf inf pi];
    StartVal = [repmat(Radius',1,2) repmat(1.5,size(Radius,2),1)]; % Use radius from find_cell as a start point
elseif Par.centering
    FitEqn = @(a, b, phi, dx, dy, x) sqrt((a * cos(0.0175*x + phi) + dx).^2 + (b * sin(0.0175*x + phi) + dy).^2 );
    lb = [0 0 0 -size(Imstack{1}{1,1})/2];
    ub = [inf inf pi size(Imstack{1}{1,1})/2];
    StartVal = [repmat(Radius',1,2) repmat([1.5 0 0],size(Radius,2),1)];
else
    error('Centering argument must have value 0 or 1')
end

% Perform unwrapping and fitting with updated centre locations
[Fits, Ia, Unwrapped] = UnwrapAndFit(Imstack, FitEqn, Radius,  Centres + Offset(2:3,:), lb, ub, StartVal, Par);

% Fix the equivalent fitting equations problem - if b>a
Fits(3,:) = Fits(3,:) + (Fits(2,:) > Fits(1,:)) * pi/2;
Fits(3,:) = mod(Fits(3,:)+pi/2, pi) - pi/2;
Fits(1:2,:) = [max(Fits(1:2,:)); min(Fits(1:2,:))];
switch nargout
    case 3; varargout = {Unwrapped, Ia};
    case 4; varargout = {Unwrapped, Ia, FitEqn};
    case 5; varargout = {Unwrapped, Ia, FitEqn, Offset};
end
end

function [Fits, Ia, Unwrapped] = UnwrapAndFit(Imstack, FitEqn, Radius, Centres, lb, ub, StartVal, Par)

    % Find the fitting variables to determine size of fits array
    FitVars = coeffnames(fittype(FitEqn));

    % Preallocate
    Theta = linspace(1,360,Par.n_theta);
    Rs = (1:round(Par.sc_up * max(Radius)))';
    Fits = zeros(length(FitVars),length(Imstack{1}));

    % For speed, everything is in one call to interp3, returning single
    % which is then cast to uint16. The first argument is the whole
    % imagestack, the remaining are x, y, z query points, all arguments are
    % cast to single for memory efficiency
    Unwrapped = uint16(interp3(...
        single(cat(3,Imstack{1}{:,1})),...                                                              % V
        single(reshape(Centres(1,:),1,1,length(Centres)) + Rs.*cos(Theta*pi/180)),...                    % Xq
        single(reshape(Centres(2,:),1,1,length(Centres)) + Rs.*sin(Theta*pi/180)),...                    % Yq
        single(repmat(reshape(1:length(Imstack{1}),1,1,length(Imstack{1})),length(Rs),Par.n_theta)),...  % Zq
        Par.inter_method));                                                                             % METHOD
    fprintf('Finished unwrapping at %gs\n',toc)

    % Find the maxes and apply filtering based on deviation from median
    % (will need large tolerance for large deformation) - maybe 0.25 for D
    % = 0.1 To make this more robust, only check for maxes away from the
    % centre (r=0) - do this by indexing a range restricted according to
    % Par.sc_down
    [~, Ia] = max(Unwrapped(floor(Par.sc_up*max(Radius)*Par.sc_down):end,:,:));
    Ia = Ia + floor(Par.sc_up*max(Radius)*Par.sc_down);
    idxa = Ia > (1-Par.tol)*median(Ia,2) & Ia < (1+Par.tol)*median(Ia,2);
    
    % For each frame, take angles corresponding to fit points and perform
    % the fit
    disp('Starting fitting')
    for frame = 1:length(Imstack{1})
        % Take corresponding angles, repeat them and unwrap
        th_fit = repmat(Theta(idxa(:,:,frame)),1,Par.n_reps) + ...
            reshape(360*(0:Par.n_reps-1).*ones(length(Theta(idxa(:,:,frame))),1),1,[]);
        fitobj = fit(double(th_fit'),repmat(Ia(:,idxa(:,:,frame),frame)',Par.n_reps,1),FitEqn,...
            'Lower', lb, 'Upper', ub, 'Start', StartVal(frame,:));
        for VarN = 1:length(FitVars) % Extract fitted values from cfit object
            Fits(VarN, frame) = fitobj.(FitVars{VarN});
        end
        prog = ceil(100 * frame / length(Imstack{1}));
        fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
    end
    fprintf('%s\r',repmat(' ',1,104))
    fprintf('Fitted data at %gs\n',toc)
end