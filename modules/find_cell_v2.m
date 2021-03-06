function [out, par] = find_cell_v2(Imstack, varargin)
%[ out ] = find_cell( Imstack ) - finds cell in image
% Use Gaussian filter, Laplacian for edge enhancement, then use Hough 
%   circles to find circles. Iterates over frames
% Returns a struct array with fields centres, radius, fails, crop.
% If multiple circles are detected, return the closest to the centre
% Fails contains 1 for each frame that has no circles, 2 for multiple
%   circles, 0 otherwise
%
% Additional option can be supplied in name-value pairs. Possibilities
%   are:
%   subFoV  -   sub Field of View to analyze - in form [y1; y2; x1; x2]
%   Gfilt   -   sigma value for Gaussian filter after edge detection
%   Lsigma  - Sigma value for Laplacian filter (edge threshold level)
%   Lalpha  - Alpha value for Laplacian filter (detail smoothing level)
%   Lbeta   - Beta value for Laplacian filter (dynamic range enhancement)
%   Rs      -   [min_R, max_R] cell radius range to check in
%   Polarity-   find a bright or dark circle
%   EThresh -   edge threshold for imfindcircles (using this causes ETScale to be ignored)
%   ETScale -   scale factor for frame-by-frame edge threshold
%   Sensitivity - sensitivity of imfindcircles ([0 1], higher finds more circles)
%   max_d   -   maximum distance from the centre to the centre of the cell
%   sc_up   -   radius multiplier for making crop box

ImW = size(Imstack{1}{1,1},2);
ImH = size(Imstack{1}{1,1},1);
if ImW >= 1920
    def_FoV = ceil([[3/16; 11/16] * ImH; [1/4; 3/4] * ImW]);
else
    def_FoV = [1; ImH; 1; ImW];
end

% Keep field names and defaults up to date here
fields = {'Gfilt',  'Lsigma', 'Lalpha', 'Lbeta', 'max_d', 'Rs', ...
    'subFoV', 'sc_up', 'Polarity','EThresh','ETScale','Sensitivity'...
    'Verbose'};
defaults = {3, 0.1, 5, 10, 250, [100 200], ...
    def_FoV, 1.25, 'dark',NaN, 1, 0.85,...
    false};

% Create parameters struct
par = cell2struct(defaults, fields,2);

% Parse inputs and create struct with parameters given
if nargin > 1
    if mod(nargin,2) == 0
        error('Please supply arguments in name-value pairs');
    end
    disp('Input arguments for find_cell:')
    for field = 1:(nargin - 1)/2
        if size(varargin{2 * field},1) ~= 0
            if size(varargin{2 * field}, 1) == 1
                disp([varargin{2 * field - 1}, ' = ', num2str(varargin{2 * field})]);
            else
                disp([varargin{2 * field - 1}, ' = ', num2str(varargin{2 * field}')]);
            end
            par.(varargin{2*field - 1}) = varargin{2*field};
        end
    end
end

tic
% Allocate output struct array
fields = {'centres', 'radius', 'find_fails', 'crop'};
defaults = {[0; 0], 0, 0, [NaN; NaN; NaN; NaN]};
out(1:size(Imstack{1},1)) = cell2struct(defaults, fields, 2);

% Allocate working cells
Ccell = cell(size(Imstack{1}, 1), 1);
Rcell = cell(size(Imstack{1}, 1), 1);

% Iterate over frames
for frame = 1:size(Imstack{1}, 1)
    prog = ceil(100 * frame / size(Imstack{1},1));
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
    % Gaussian filter, Laplacian filter, Hough circles 
    % Use indexing to pull out the sub field of view - the cell is never
    % near the edge.

    ImFilt = imgaussfilt(Imstack{1}{frame,1} ...
        (par.subFoV(1):par.subFoV(2), par.subFoV(3):par.subFoV(4))...
        , par.Gfilt);
    ImLap = locallapfilt(ImFilt,par.Lsigma, par.Lalpha, par.Lbeta);
    
    % Find circular things - if user hasn't supplied edge threshold value,
    % use the default from graythresh.
    if ~isfinite(par.EThresh)
        EThresh = graythresh(imadjust(ImLap)) * par.ETScale;
    else
        EThresh = par.EThresh;
    end
    [Cs, Rs] = imfindcircles(imadjust(ImLap), par.Rs,  'ObjectPolarity',par.Polarity,'EdgeThreshold',EThresh,'Sensitivity',par.Sensitivity);
    
    % We are only interested in circles that are close to the centre
    % Go through the found circles, and add the good ones to a "keep" array
    keepC = [];
    keepR = [];
    for circ = 1:size(Rs, 1)
        % For each circle, only keep it if it's less than max_d from centre
        % Offset by the subFoV
        if sqrt( sum( (Cs(circ,:) - ...
                [ImW/2 - par.subFoV(3), ImH/2 - par.subFoV(1)]).^2)) ...
                < par.max_d
            keepC = [keepC; Cs(circ, :)];
            keepR = [keepR; Rs(circ)];
        end
    end
    fprintf('%s\r',repmat(' ',1,102))
    if par.Verbose
        if size(keepR,1) == 0
            disp(['No circles detected in frame ', num2str(frame)])
            if size(Rs,1) > 0
                fprintf('Discarded %i circles due to max distance from centre > %i\n',...
                    size(Rs,1), par.max_d)
            end
        elseif size(keepR,1) == 1
            disp(['One circle detected in frame ', num2str(frame)])
        else
            disp([num2str(size(keepR,1)), ' circles detected in frame ', ...
                num2str(frame)])
        end
    end
    Ccell{frame} = keepC;
    Rcell{frame} = keepR;
end
 
n_circles = zeros(size(Imstack{1},1), 1);

for frame = 1:size(Imstack{1},1)
    n_circles(frame) = size(Rcell{frame},1);
end

allC = zeros(sum(n_circles), 2);
allR = zeros(sum(n_circles), 1);

for frame = 1:size(Imstack{1},1)
    if frame == 1
        idx = 1;
    else
        idx = sum(n_circles(1:frame-1)) + 1;
    end
    allC(idx:idx - 1 + n_circles(frame), :) = Ccell{frame};
    allR(idx:idx -1 + n_circles(frame)) = Rcell{frame};
end

% The cell's position doesn't change much over the course of the imaging.
% Taking the mode should make it possible to find the correct cell.
% Round to nearest 10px to account for noise.
modeC = mode(round(allC, -1));

% After finding circles, we want the circles closest to the modal position
for frame = 1:size(Imstack{1},1)
    
    % If multiple circles are found in this frame
    if size(Rcell{frame}, 1) > 1
        % Find the best circle by choosing the one closest to the mode
        % bsxfun subtracts the mode from each row, which is then squared
        % and summed. The minimum value (distance to mode) is kept.
        [~, I] = min(sum(bsxfun(@minus, Ccell{frame}, modeC).^2,2));
        
        % Return the best circle
        out(frame).centres = Ccell{frame}(I,:)';
        out(frame).radius = Rcell{frame}(I);
        out(frame).find_fails = 2;
        
    % If no circles are detected: return NaN
    elseif size(Rcell{frame}, 1) < 1
        out(frame).centres(:) = NaN;
        out(frame).radius = NaN;
        out(frame).find_fails = 1;
        
    % Exactly 1 circle detected: return that circle
    else
        out(frame).centres = Ccell{frame}';
        out(frame).radius = Rcell{frame};
        out(frame).find_fails = 0;
    end
    
    % Add the offset from cropping
    out(frame).centres = out(frame).centres + ...
        [par.subFoV(3); par.subFoV(1)]; 
    % Create the crop values for this cell
    if out(frame).find_fails ~= 1
        out(frame).crop(1:2) = out(frame).centres(1) + ...
            [-1; 1] * par.sc_up * out(frame).radius;
        out(frame).crop(3:4) = out(frame).centres(2) + ...
            [-1; 1] * par.sc_up * out(frame).radius;
        for idx = 1:4
            if out(frame).crop(idx) < 1
                out(frame).crop(idx) = 1;
            end
        end
    else
        out(frame).crop = [NaN; NaN; NaN; NaN];
    end
    % Crop is used for indexing later, so must be integer valued
    out(frame).crop = round(out(frame).crop);
end

toc
end

