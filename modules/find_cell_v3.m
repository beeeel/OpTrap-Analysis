function [out] = find_cell_v3(Im, ImW, ImH, varargin)
%[ out ] = find_cell_v3( Im, ImW, ImH ) - finds cell in image
% Use Im = cat(2,Imstack{1}{:,1}) to convert from OME cell array.
%
% Use Gaussian filter, Laplacian for edge enhancement, then use Hough 
%   circles to find circles. Vectorized circle detection, iterative post
%   selection.
% Returns a struct array with fields centres, radius, fails, crop.
% If multiple circles are detected, return the closest to the centre
% Fails contains 1 for each frame that has no circles, 2 for multiple
%   circles, 0 otherwise
%
% Additional option can be supplied in name-value pairs. Possibilities
%   are:
%   Gfilt   -   sigma value for Gaussian filter after edge detection
%   Lsigma  -   Sigma value for Laplacian filter (edge threshold level)
%   Lalpha  -   Alpha value for Laplacian filter (detail smoothing level)
%   Lbeta   -   Beta value for Laplacian filter (dynamic range enhancement)
%   Rs      -   [min_R, max_R] cell radius range to check in (integer values)
%   Polarity-   find a bright or dark circle
%   max_d   -   maximum distance from the centre to the centre of the cell
%   sc_up   -   radius multiplier for making crop box
%   weight  -   Weighting for circle selection (0 - 1)
%
% The imfindcircles call can be sped up by pre-measuring the radius
% (manually) and giving good values for Rs (radius range to search within)
%
% weight is used to balance between the quality metric from imfindcircles,
% and the modal position based metric. Around 0.01-0.1 seems appropriate
% (this is a total guess, I might come and update this comment at the end
% of my PhD)

N_frames = size(Im, 2)/ImW;
if N_frames ~= floor(N_frames)
    error('Value of ImW %i does not match size of input image %i',ImW, size(Im,2))
end

% Keep field names and defaults up to date here
fields = {'Gfilt',  'Lsigma', 'Lalpha', 'Lbeta', 'max_d', 'Rs', ...
    'sc_up', 'Polarity'};
defaults = {10, 0.05, 0.1, 100, 250, [100 200], ...
    1.25, 'dark'};

% Create parameters struct
par = cell2struct(defaults, fields,2);

% Parse inputs and create struct with parameters given
if nargin > 3
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
out(1:N_frames) = cell2struct(defaults, fields, 2);

% Allocate working cells
Ccell = cell(N_frames, 2);
Rcell = cell(N_frames, 1);
Mcell = cell(N_frames, 1);

% Preprocess all frames - concatenate columns for fast indexing while
% having a 2D plane for image functions

% Gaussian filter, Laplacian filter, Hough circles
% Use indexing to pull out the sub field of view - the cell is never
% near the edge.
disp('Processing images')
Im = imgaussfilt(Im, par.Gfilt);
Im = locallapfilt(Im,par.Lsigma, par.Lalpha, par.Lbeta);
disp('Finding circles')
[Cs, Rs, Ms] = imfindcircles(imadjust(Im), par.Rs,  'ObjectPolarity',par.Polarity);
disp('Post-selecting circles')
disp('Applying distance constraint')
% We are only interested in circles that are close to the centre
% Go through the found circles, and add the good ones to a "keep" array
for circ = 1:size(Rs, 1)
    prog = ceil(100 * circ / size(Rs,1));
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
    % For each circle, only keep it if it's less than max_d from centre
    % Wrap centres to the frame you're looking in using floor.
    if sqrt( sum( ([floor(Cs(circ,1)/ImW)*ImW, Cs(circ,2)] - ...
            [ImW/2, ImH/2]).^2 ) ) ...
            < par.max_d
        % Calculate which frame this circle is centred within and place it
        % in cell array
        frame = floor(Cs(circ,1)/ImW) + 1;
        Ccell{frame}(end+1,:) = [mod(Cs(circ, 1), ImW), Cs(circ,2)];
        Rcell{frame}(end+1) = Rs(circ);
        Mcell{frame}(end+1) = Ms(circ);
    end
end
fprintf('%s\r',repmat(' ',1,106))

disp('Selecting best circles for output')
% The cell's position doesn't change much over the course of the imaging.
% Taking the mode should make it possible to find the correct cell.
% Round to nearest 10px to account for noise.
modeC = mode(round(Cs, -1));

% After finding circles, we want the circles closest to the modal position
for frame = 1:N_frames
    prog = ceil(100 * frame / N_frames);
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
    % If multiple circles are found in this frame
    if size(Rcell{frame}, 1) > 1
        % Find the best circle by choosing the one closest to the mode
        % bsxfun subtracts the mode from each row, which is then squared
        % and summed. This is added to the metric from imfindcircles, with
        % weighting defined by weight. The minimum value is kept.
        [~, I] = min(weight * sum(bsxfun(@minus, Ccell{frame}, modeC).^2,2) ...
            + (1-weight) * Mcell{frame});
        
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
    out(frame).centres = out(frame).centres; % Create the crop values for this cell
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
fprintf('%s\r',repmat(' ',1,106))

end