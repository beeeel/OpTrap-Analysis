%% Measure cell by radial unwrapping
% Start with the centre of the cell, as found by find_cell, then draw lines
% radially to create line plots of intensity (interpolating as necessary).
% Find edge of cell from intensity plots, and then fit to parametric
% equation for ellipse.
%% Load
% Load imstack and info file - this contains the results from find_cell
CellType = 'LS174T';
Set = 'normoxia';
Num = '11';
[Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num);
%% Unpack stuff from info
% These were made by find_cell_v2
centres = [info.mCentres];
radii = min(abs([info.mCentres] - min(size(Imstack{1}{1,1}))));
%% Unwrap and fit to radial equation for ellipse
% Parameters

fields = {'sc_up', 'n_theta', 'n_reps', 'tol', 'inter_method', 'sc_down',...
    'centering', 'ifNaN','parallel','edge_method'};
defaults = {1, 360, 5, 0.15, 'linear', 0.5,...
    0, 'mean',false, 'gradient'};
Par = cell2struct(defaults, fields,2);

% Parametric eqn of off-centre circle (x = cos(t) + dx, y = sin(t) + dy),
% convert to polar (r = x.^2 + y.^2)
CircleEqn = @(r, dx, dy, x) sqrt((r .* cos(0.0175.*x) + dx).^2 + (r .* sin(0.0175.*x) + dy).^2);
lb = [0, -size(Imstack{1}{1,1})/2];
ub = [Par.sc_up * max(radii), size(Imstack{1}{1,1})/2];
StartVal = [radii', repmat([0, 0],size(radii,2), 1)];

% Perform unwrapping and fitting for off-centreness
[Offset, ~, ~, ~] = UnwrapAndFit(Imstack, CircleEqn, radii,  centres, lb, ub, StartVal, Par);

% From parametric equation of ellipse (x = a cos(t), y = b sin(t)), take to
% polar, r = (a.cos(f * theta + phi))^2 +  (b.sin(f * theta + phi))^2,
% where f is the frequency - twice per 2pi, theta is the angle (x below),
% and phi is the phase offset (orientation of ellipse)
fiteqn = @(a, b, c, x) (a.^2) *cos(0.0175*x + c).^2 + (b.^2)*sin(0.0175*x + c).^2 ;
lb = [0 0 0];
ub = [inf inf pi];
StartVal = [repmat(radii',1,2) repmat(1.5,size(radii,2),1)]; % Use radius from find_cell as a start point

tic
% For memory efficiency, everything is in a single call to interp3, which
% is then cast to uint16. The first argument is the whole imagestack, the
% remaining are x, y, z query points, all arguments are cast to single for
% memory efficiency
[fits, Ia, Unwrapped, Errs] = UnwrapAndFit(Imstack, fiteqn, radii, centres, lb, ub, StartVal, Par);
% Fix the equivalent fitting equations problem - if b>a
fits(3,:) = fits(3,:) + (fits(2,:) > fits(1,:)) * pi/2;
fits(3,:) = mod(fits(3,:)+pi/2, pi) - pi/2;
fits(1:2,:) = [max(fits(1:2,:)); min(fits(1:2,:))];
fprintf('%s\n',repmat(' ',1,104))
fprintf('Fitted data in %gs\n',toc)
%% Get gradient image - it's a race!
tic
GradStack = Imstack;
for fr = 1:length(GradStack{1})
    [Gmag, ~] = imgradient(Imstack{1}{fr,1});
    GradStack{1}{fr,1} = Gmag;
end
toc
%% The second contender: vectorisation!
tic
Kernel = fspecial('sobel');
Ims = single(cat(3,Imstack{1}{:,1}));
Gy = imfilter(Ims, -Kernel);
Gx = imfilter(Ims, -Kernel');
Gmax = Gx.^2 + Gy.^2;
toc
%%
[Gfits, ~, GUnwrapped, ~] = UnwrapAndFit(GradStack, fiteqn, radii, centres, lb, ub, StartVal, Par);

%% Unwrap functions
function [Fits, Ia, Unwrapped, Errs] = UnwrapAndFit(Imstack, FitEqn, Radius, Centres, lb, ub, StartVal, Par)

    % Find the fitting variables to determine size of fits array
    FitVars = coeffnames(fittype(FitEqn));

    % Preallocate
    Theta = linspace(1,360,Par.n_theta);
    Rs = (1:round(Par.sc_up * max(Radius)))';
    Fits = zeros(length(FitVars),length(Imstack{1}));
    Errs = Fits;

    % Debugging memory usage
    Sz = size(Imstack{1}{1,1});
    Frs = length(Imstack{1});
    Sgl = single(1);
    Sgl = whos('Sgl');
    SglMem = Sgl.bytes;
    MemNeeded = (Sz(1) * Sz(2) * Frs + 4 * ( Par.n_theta * length(Rs) * Frs)) * SglMem ./ 1e9;
    fprintf('About to request %g GB of memory. Gulp!\n',MemNeeded)
    
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
    [Ia, idxa] = FindMaxes(Unwrapped, Imstack, Radius, Par);
    
    % For each frame, take angles corresponding to fit points and perform
    % the fit
    disp('Starting fitting')
    if ~Par.parallel
        for frame = 1:length(Imstack{1})
            % Take corresponding angles, repeat them and unwrap
            th_fit = repmat(Theta(idxa(:,:,frame)),1,Par.n_reps) + ...
                reshape(360*(0:Par.n_reps-1).*ones(length(Theta(idxa(:,:,frame))),1),1,[]);
            fitobj = fit(double(th_fit'),repmat(Ia(:,idxa(:,:,frame),frame)',Par.n_reps,1),FitEqn,...
                'Lower', lb, 'Upper', ub, 'Start', StartVal(frame,:));
            for VarN = 1:length(FitVars) % Extract fitted values from cfit object
                Fits(VarN, frame) = fitobj.(FitVars{VarN});
            end
            ConfInt = confint(fitobj);
            Errs(:,frame) = (ConfInt(2,:) - ConfInt(1,:))/2;
            prog = ceil(100 * frame / length(Imstack{1}));
            fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
        end
    else
        
        n_reps = Par.n_reps;
        n_frames = length(Imstack{1});
        
        Ias = cell(n_frames,1);
        Thetas = cell(n_frames,1);
        for frame = 1:n_frames
            Thetas{frame} = Theta(idxa(:,:,frame));
            Ias{frame} = Ia(:,idxa(:,:,frame))';
        end
        
        parfor frame = 1:length(Imstack{1})
            % Take corresponding angles, repeat them and unwrap
            th_fit = repmat(Thetas{frame},1,n_reps) + ...
                reshape(360*(0:n_reps-1).*ones(length(Thetas{frame}),1),1,[]);
            fitobj = fit(double(th_fit'),repmat(Ias{frame},n_reps,1),FitEqn,...
                'Lower', lb, 'Upper', ub, 'Start', StartVal(frame,:));
            
            Fits(:,frame) = coeffvalues(fitobj);
            ConfInt = confint(fitobj);
            Errs(:,frame) = (ConfInt(2,:) - ConfInt(1,:))/2;
        end
    end
    fprintf('%s\r',repmat(' ',1,104))
    fprintf('Fitted data at %gs\n',toc)
end

function [Ia, idxa] = FindMaxes(Unwrapped, Imstack, Radius, Par)
switch Par.edge_method
    case 'simple'
        [~, Ia] = max(Unwrapped(floor(Par.sc_up*max(Radius)*Par.sc_down):end,:,:));
        Ia = Ia + floor(Par.sc_up*max(Radius)*Par.sc_down);
        idxa = Ia > (1-Par.tol)*median(Ia,2) & Ia < (1+Par.tol)*median(Ia,2);
    case 'RadialGradient'
        Unwrapped = Unwrapped(floor(Par.sc_up*max(Radius)*Par.sc_down):end,:,:);
        Gy = zeros(size(Unwrapped));
        Filt = zeros(size(Unwrapped));
        tic
        for frame = 1:size(Unwrapped,3)
            Filt(:,:,frame) = imgaussfilt(Unwrapped(:,:,frame),5,'FilterSize',[31,1]);
            [~, Gy(:,:,frame )] = imgradientxy(Unwrapped(:,:,frame));
        end
        toc
    case 'MaxGradient'
       
        
end
end
