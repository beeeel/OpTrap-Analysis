function [Fits, varargout] = unwrap_cell_v4(Imstack, Centres, Radius, varargin)
%unwrap_cell_v4 - perform radial unwrapping by interpolation and fit to 2D
%ellipse equation. Handle off-centre by fitting twice
% [Fits, varargout] = unwrap_cell_v4(Imstack, Centres, Radius, varargin) -
% if Imstack has N frames, centres and radii must be 2xN and 1xN and
% contain [X_centre, Y_centre] and radius in each column, respectively.
%
% Varargout = {Unwrapped, Ia, FitEqn, Offset, FitErrs}
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
    'centering', 'ifNaN','parallel','weighted'};
defaults = {1.2, 360, 5, 0.15, 'linear', 0.5,...
    0, 'mean',false, true};

tic

Par = ParseInputs(fields, defaults, varargin{:});

% Fix any NaNs from find_cell failing
[Centres, Radius] = FixNaNs(Centres, Radius, Imstack, Par);

% Get the circle equation
[CircleEqn, lb, ub, StartVal] = GetEqn('circle', Imstack, Par, Radius);

% Debugging memory usage
Var = whos;
fprintf('Currently hogging %g GB of memory! Whoops!\n',sum([Var.bytes])./1e9)
clear Var;

% Perform unwrapping and fitting for off-centreness
[Offset, ~, ~] = UnwrapAndFit(Imstack, CircleEqn, Radius,  Centres, lb, ub, StartVal, Par);

% Get the ellipse equation
[FitEqn, lb, ub, StartVal] = GetEqn('ellipse', Imstack, Par, Radius);

% Perform unwrapping and fitting with updated centre locations
[Fits, Unwrapped, FitErrs] = UnwrapAndFit(Imstack, FitEqn, Offset(1,:),  Centres + Offset(4:5,:), lb, ub, StartVal, Par);

% Fix the equivalent fitting equations problem - if b>a
Fits(3,:) = Fits(3,:) + (Fits(2,:) > Fits(1,:)) * pi/2;
Fits(3,:) = mod(Fits(3,:)+pi/2, pi) - pi/2;
Fits(1:2,:) = [max(Fits(1:2,:)); min(Fits(1:2,:))];
switch nargout
    case 2; varargout = {Unwrapped};
    case 3; varargout = {Unwrapped, FitEqn};
    case 4; varargout = {Unwrapped, FitEqn, Offset};
    case 5; varargout = {Unwrapped, FitEqn, Offset, FitErrs};
    case 6; varargout = {Unwrapped, FitEqn, Offset, FitErrs, Par};
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Functions start here
    function [Fits, Unwrapped, Errs] = UnwrapAndFit(Imstack, FitEqn, Radius, Centres, lb, ub, StartVal, Par)
        %% Fairly self-explanatory tbh. Unwraps data using interp3, finds edges, and fits to fit equation
        % Find the fitting variables to determine size of fits array
        FitVars = coeffnames(fittype(FitEqn,'independent',{'Theta','R'}));
        N_Frs = length(Imstack{1});
        
        % Preallocate
        Theta = linspace(1,360,Par.n_theta);
        Rs = (1:round(Par.sc_up * max(Radius)))';
        Fits = zeros(length(FitVars),length(Imstack{1}));
        Errs = Fits;
        
        % Debugging memory usage
        Sz = size(Imstack{1}{1,1});
        Sgl = single(1);
        Sgl = whos('Sgl');
        SglMem = Sgl.bytes;
        MemNeeded = (Sz(1) * Sz(2) * N_Frs + 4 * ( Par.n_theta * length(Rs) * N_Frs)) * SglMem ./ 1e9;
        fprintf('About to request %g GB of memory. Gulp!\n',MemNeeded)
        
        % For speed, everything is in one call to interp3, returning single
        % which is then cast to uint16. The first argument is the whole
        % imagestack, the remaining are x, y, z query points, all arguments are
        % cast to single for memory efficiency
        Unwrapped = double(interp3(...
            single(cat(3,Imstack{1}{:,1})),...                                                              % V
            single(reshape(Centres(1,:),1,1,N_Frs) + Rs.*cos(Theta*pi/180)),...                    % Xq
            single(reshape(Centres(2,:),1,1,N_Frs) + Rs.*sin(Theta*pi/180)),...                    % Yq
            single(repmat(reshape(1:N_Frs,1,1,N_Frs),length(Rs),Par.n_theta)),...  % Zq
            Par.inter_method));                                                                             % METHOD
        fprintf('Finished unwrapping at %gs\n',toc)
        
        [Fits, Errs] = N_DoUnwrappedFits(Theta, Rs, Unwrapped, FitEqn, lb, ub, StartVal, Fits, Errs, N_Frs, Par);
        
        fprintf('%s\r',repmat(' ',1,104))
        fprintf('Fitted data at %gs\n',toc)
    end

    function [Fits, Errs] = N_DoUnwrappedFits(Theta, Rs, Unwrapped, FitEqn, lb, ub, StartVal, Fits, Errs, N_Frs, Par)
        %% Performs fits - parallel or non, or weighted
        % For each frame, take angles corresponding to fit points and perform
        % the fit
        disp('Starting fitting')
        ThisFit = fittype(FitEqn,'independent',{'Theta','R'});
        FitVars = coeffnames(ThisFit);
        if ~Par.parallel
            if Par.weighted                
                WFun = @(r, theta, frame) repmat(normalize(normpdf(1:max(Rs),Radius(frame), Radius(frame)/2),'scale','first')',1,length(theta));
                for frame = 1:N_Frs
                    % Take corresponding angles, repeat them and unwrap
                    Weights = WFun(Rs, Theta, frame);
                    [ThOut, rOut, IOut, WOut] = prepareSurfaceData(Theta, Rs, Unwrapped(:,:,frame), Weights);
                    fitobj = fit([ThOut, rOut], IOut,ThisFit,'Upper',ub,'Lower',lb,'Start',StartVal(frame,:),'Weights',WOut);
                    for VarN = 1:length(FitVars) % Extract fitted values from cfit object
                        Fits(VarN, frame) = fitobj.(FitVars{VarN});
                    end
                    ConfInt = confint(fitobj);
                    Errs(:,frame) = (ConfInt(2,:) - ConfInt(1,:))/2;
                    ProgressBar(frame./N_Frs)
                end
            else
                for frame = 1:N_Frs
                    % Take corresponding angles, repeat them and unwrap
                    [ThOut, rOut, IOut] = prepareSurfaceData(Theta, Rs, Unwrapped(:,:,frame));
                    fitobj = fit([ThOut, rOut], IOut,ThisFit,'Upper',ub,'Lower',lb,'Start',StartVal(frame,:));
                    for VarN = 1:length(FitVars) % Extract fitted values from cfit object
                        Fits(VarN, frame) = fitobj.(FitVars{VarN});
                    end
                    ConfInt = confint(fitobj);
                    Errs(:,frame) = (ConfInt(2,:) - ConfInt(1,:))/2;
                    ProgressBar(frame./N_Frs)
                end
            end
        else
            error('Parallel not programmed yet')
            %{
    if Par.weighted
        warning('Weighted argument ignored because running parallel. Write some more code!')
    end
    n_reps = Par.n_reps;
    
    Ias = cell(N_Frs,1);
    Thetas = cell(N_Frs,1);
    for frame = 1:N_Frs
        Thetas{frame} = Theta(idxa(:,:,frame));
        Ias{frame} = Ia(:,idxa(:,:,frame))';
    end
    
    parfor frame = 1:N_Frs
        % Take corresponding angles, repeat them and unwrap
        th_fit = repmat(Thetas{frame},1,n_reps) + ...
            reshape(360*(0:n_reps-1).*ones(length(Thetas{frame}),1),1,[]);
        fitobj = fit(double(th_fit'),repmat(Ias{frame},n_reps,1),FitEqn,...
            'Lower', lb, 'Upper', ub, 'Start', StartVal(frame,:));
        
        Fits(:,frame) = coeffvalues(fitobj);
        ConfInt = confint(fitobj);
        Errs(:,frame) = (ConfInt(2,:) - ConfInt(1,:))/2;
    end
            %}
        end
    end
    function [FitEqn, lb, ub, StartVal] = GetEqn(Func, Imstack, Par, Radius)
        %% Get an equation and bounds for fitting
        switch Func
            case 'circle'
                % Parametric eqn of off-centre circle (x = cos(t) + dx, y = sin(t) + dy),
                % convert to polar (r = x.^2 + y.^2)
                FitEqn =@(r, sigma, amp, dx, dy, Theta, R) ...
                    amp .* (R - sqrt((r .* cos(0.0175.*Theta) + dx).^2 + (r .* sin(0.0175.*Theta) + dy).^2)) .* ...
                    exp(-((R - sqrt((r .* cos(0.0175.*Theta) + dx).^2 + (r .* sin(0.0175.*Theta) + dy).^2)).^2)./(2 .* sigma .^2)) ./ ...
                    (sqrt(2 .* pi) * sigma.^3);
                lb = [0, 0, 0, -size(Imstack{1}{1,1})/2];
                ub = [Par.sc_up * max(Radius), inf, inf, size(Imstack{1}{1,1})/2];
                StartVal = [Radius', repmat([5, 1e4, 0, 0],size(Radius,2), 1)];
            case 'ellipse'
                % From parametric equation of ellipse (x = a cos(t), y = b sin(t)), take to
                % polar, r = (a.cos(f * theta + phi))^2 +  (b.sin(f * theta + phi))^2,
                % where f is the frequency - once per 2pi (doubled by squaring), theta is
                % the angle (x below is in degrees), and phi is the phase offset
                % (orientation of ellipse). For an off-centre elipse, x = a cos(t) + dx, y
                % = a sin(t) + dy.
                
                % Take the parametric equation and put it into the derivative of
                % gaussian formula to create a surface equation for a centred
                % ellipse
                if ~Par.centering
                    FitEqn = @(a, b, phi, sigma, amp, Theta, R) ...
                        amp .* (R - sqrt((a .* cos(0.0175.*Theta + phi)).^2 + (b .* sin(0.0175.*Theta + phi)).^2)) .* ...
                        exp(-((R - sqrt((a .* cos(0.0175.*Theta + phi)).^2 + (b .* sin(0.0175.*Theta + phi)).^2)).^2)./(2 .* sigma .^2)) ./ ...
                        (sqrt(2 .* pi) * sigma.^3);
                    lb = [0, 0, 0, 0, 0];
                    ub = [inf, inf, pi, 100, inf];
                    StartVal = [repmat(Radius',1,2) repmat([1.5, 5, 1e4],size(Radius,2),1)]; % Use radius from find_cell as a start point
                elseif Par.centering
                    FitEqn =  @(a, b, phi, sigma, amp, dx, dy, Theta, R) ...
                        amp .* (R - sqrt((a .* cos(0.0175.*Theta + phi) + dx).^2 + (b .* sin(0.0175.*Theta + phi) + dy).^2)) .* ...
                        exp(-((R - sqrt((a .* cos(0.0175.*Theta + phi) + dx).^2 + (b .* sin(0.0175.*Theta + phi) + dy).^2)).^2)./(2 .* sigma .^2)) ./ ...
                        (sqrt(2 .* pi) * sigma.^3);
                    lb = [0, 0, 0, 0, 0, -size(Imstack{1}{1,1})/2];
                    ub = [inf, inf, pi, inf, inf, size(Imstack{1}{1,1})/2];
                    StartVal = [repmat(Radius',1,2) repmat([1.5 5 1e4 0 0],size(Radius,2),1)];
                else
                    error('Centering argument must be true or false')
                end
            otherwise
                error('miscall to GetEqn - you need to specify which equation you want')
        end
    end

    function [Centres, Radius] = FixNaNs(Centres, Radius, Imstack, Par)
        %% Fix NaN values in centres and radius arrays
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
                % Assume the cell is centred and fills  the FoV
                [ImH, ImW] = size(Imstack{1}{1,1});
                Centres(:,~isfinite(Radius)) = repmat([ImH/2 ImW/2],1,sum(~isfinite(Radius)));
                Radius(~isfinite(Radius)) = (ImW < ImH) * ImW/2 + (ImW > ImH) * ImH/2;
            otherwise
                error('ifNaN must have value ''mean'' or ''last'' or ''centre''')
        end
    end

    function Par = ParseInputs(fields, defaults, varargin)
        %% Parse inputs
        Par = cell2struct(defaults, fields,2);
        def_argin = 2;
        
        % Parse inputs and create struct with parameters given
        if ~isempty(varargin)
            if mod(nargin,2) == 1 % Imstack counts towards nargin
                error('Please supply arguments in name-value pairs');
            end
            stack = dbstack;
            % For each field, display it depending on its size and contents
            fprintf('Input arguments for %s:\n',stack(2).name)
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
    end
end