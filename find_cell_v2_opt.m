function [ met ] = find_cell_v2_opt( Imstack , GfiltIn, LsigmaIn, LalphaIn, LbetaIn)
%[ met ] = find_cell_v2_opt( Imstack ) - measures parameter optimisation
% Use large filter Laplacian of Guassian (log) for edge detection, filter
%   with gaussian filter, then use Hough circles to find circles.
% Returns the total number of frames where at least 1 good circle was found

% Keep field names and defaults up to date here
fields = {'Gfilt',  'Lsigma', 'Lalpha', 'Lbeta', 'max_d', 'Rs', ...
    'subFoV', 'sc_up', 'Polarity'};
defaults = {10, 250, 0.9, 2, 250, [100 200], ...
    [282; 629; 759; 1229], 1.1, 'dark'};

% Create parameters struct
par = cell2struct(defaults, fields,2);

% Take input arguments and allocate output variable
par.Gfilt = GfiltIn;
par.Lsigma = LsigmaIn;
par.Lbeta = LbetaIn;
par.Lalpha = LalphaIn;
met = 0;

tic
% Iterate over frames
for frame = 1:size(Imstack{1}, 1)
    % Edge filter, gaussian filter, Hough circles
    % Use indexing to pull out the middle quarter of the area - the cell is
    % never near the edge    
    ImFilt = imgaussfilt(Imstack{1}{frame,1} ...
        (par.subFoV(1):par.subFoV(2), par.subFoV(3):par.subFoV(4))...
        , par.Gfilt);
    ImLap = locallapfilt(ImFilt,  par.Lsigma, par.Lalpha, par.Lbeta);
    %ImBin = imbinarize(ImLap, 'global');
    
    % Find circular things
    [Cs, Rs] = imfindcircles(ImLap, par.Rs,  'ObjectPolarity',par.Polarity);
    
    % We are only interested in circles that are close to the centre
    % Go through the found circles, and add the good ones to a "keep" array
    keepN = 0;
    for circ = 1:size(Rs, 1)
        % For each circle, only keep it if it's less than max_d from centre
        % Offset by the subFoV
        if sqrt( sum( (Cs(circ,:) - ...
                [960 - par.subFoV(1,1), 540 - par.subFoV(2,1)]).^2)) ...
                < par.max_d
            keepN = keepN + 1;
        end
    end
    if keepN == 0
        met = met + 1;
    end
end


toc
end