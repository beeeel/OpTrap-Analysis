function [ met ] = find_cell_opt( Imstack , EsigmaIn, GfiltIn)
%[ met ] = find_cell( Imstack ) - measures parameter optimisation
% Use large filter Laplacian of Guassian (log) for edge detection, filter
%   with gaussian filter, then use Hough circles to find circles.
% Returns the total number of frames where at least 1 good circle was found

ImW = size(Imstack{1}{1,1},2);
ImH = size(Imstack{1}{1,1},1);

% Keep field names and defaults up to date here
fields = {'Esigma', 'Gfilt', 'edge', 'max_d', 'Rs', ...
    'subFoV', 'sc_up'};
defaults = {10, 10, 'canny', 250, [100 200], ...
    ceil([[4/16, 12/16] * ImH; [4/16, 12/16] * ImW]), 1.1};

% Create parameters struct
par = cell2struct(defaults, fields,2);

% Take input arguments and allocate output variable
par.Esigma = EsigmaIn;
par.Gfilt = GfiltIn;
met = 0;

tic
% Iterate over frames
for frame = 1:size(Imstack{1}, 1)
    % Edge filter, gaussian filter, Hough circles
    % Use indexing to pull out the middle quarter of the area - the cell is
    % never near the edge
    ImEdge = edge(Imstack{1}{frame,1} ...
        (par.subFoV(1,1):par.subFoV(1,2), par.subFoV(2,1):par.subFoV(2,2))...
        , par.edge, [], par.Esigma);
    ImFilt = imgaussfilt(double(ImEdge), par.Gfilt) * 10;
    
    % Radius should be in this range for HeLa cells in my images
    [Cs, Rs] = imfindcircles(ImFilt, par.Rs);
    
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