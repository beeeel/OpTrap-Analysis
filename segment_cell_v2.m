function [ mask ] = segment_cell_v2( Imstack, varargin)
%Segment a cell from the centre of the image
% Use computer vision to segment a cell from the image given. If no values
%  for crop are given, the cell is assumed to be centrally in the image.
% First the image is cropped, then edge detected using Canny's method,
%  Gaussian filtered and thresholded. Small objects are removed and a
%  flattening kernel is applied before a final thresholding
% Mask is returned as an array of size [image x, image y, num frames]. For
%  frames where fails(frame) == 1, mask(1,1,frame) == NaN.
% Parameters can be input as name-value pairs. Possibilities are:
%   crop    - Field of view containing the cell
%   sigma   - Standard deviation of Gaussian filter used on edge map
%   Ethresh - Threshold applied to filtered edge map
%   minSize - Minimum area for the cell (in px)
%   windowSize - Size of convolution kernel used to flatten image
%   Fthresh - Threshold applied to flattened image (final processing step)
%   ellipseFitVal - Method for mask creation
%   fails   - Vector equal in length to number of frames, containing 1 for
%               each frame to skip segmentation

% Keep field names and defaults up to date here
fields = {'crop', 'sigma', 'Ethresh', 'minSize', 'windowSize', ...
    'Fthresh', 'ellipseFitVal', 'fails'};
defaults = {[340; 740; 760; 1160], 2, 0.5, 2000, 51, ...
    0.4, 1, zeros(size(Imstack{1},1),1)};

% Parse inputs and create struct with parameters given
if nargin > 1
    if mod(nargin,2) == 0
        error('Please supply arguments in name-value pairs');
    end
    for field = 1:(nargin - 1)/2
        if size(varargin{2 * field}, 2) < 4
            disp([varargin{2 * field - 1}, ' = ', varargin{2 * field}]);
        else
            disp([varargin{2 * field - 1}, ' = size[', ...
                size(varargin{2 * field},1), size(varargin{2 * field},2)...
                , ']']);
        end
        par.(varargin{2*field - 1}) = varargin{2*field};
    end
end

% For each name in fields, ensure par has that field
for f = 1:size(fields,2)
    if ~isfield(par, fields{f})
        par.(fields{f}) = defaults{f};
    end
end

kernel = ones(par.windowSize) / par.windowSize ^ 2;

mask = zeros([ size(Imstack{1}{1,1}), size(Imstack{1},1)]);

for frame = 1:size(Imstack{1},1)
    if par.fails(frame) ~= 1
        % Crop out relevant area - if the input is a full stack, take the
        % relevant crop values from the crop array. Otherwise crop is a vector
        if size(Imstack{1},1) > 1
            subIm = Imstack{1}{frame,1}(par.crop(3,frame):par.crop(4, frame), par.crop(1, frame):par.crop(2, frame));
        else
            subIm = Imstack{1}{frame,1}(par.crop(3):par.crop(4), par.crop(1):par.crop(2));
        end
        % Edge detection using canny method
        ImEdge = edge(subIm, 'canny');
        % Apply gaussian blur and threshold to create binary map
        ImThresh = imgaussfilt(double(ImEdge), par.sigma) * 10 > par.Ethresh; % sigma = 6 was ok
        % Remove objects smaller than a certain size of pixels
        ImOpen = bwareaopen(ImThresh, par.minSize);
        % Apply a flattening/average filter and threshold
        smoothed = conv2(single(ImOpen), kernel, 'same') > par.Fthresh;
        
        
        % This is only here because it was in v1 - I only use ellipseFitVal
        % = 0. 
        if par.ellipseFitVal == 1
            %Now generate Elipse paramaters and fit using Elipse fitting function
            [semiY,semiX,centreX,centreY] = ellipse_fit_fun_final(smoothed);
            
            mask(:,:,frame) = ellipse_mask(semiX,semiY,centreX,centreY,smoothed);
            
        elseif par.ellipseFitVal ==2
            %Else if ellipseFitVal ==2, use ellipseDetection method
            % I don't think this actually works
            [bestFits] = ellipseDetection(ImEdge);
            
            mask(par.crop(3,frame):par.crop(4, frame), ...
                par.crop(1, frame):par.crop(2, frame),frame)...
                 = ellipse_mask(majorAxis,minorAxis,centreX,centreY,bwSmoothed);
        else
            %If not fitting Ellipse, just set to previous Mask (BWSmoothed)
            mask(par.crop(3,frame):par.crop(4, frame), ...
                par.crop(1, frame):par.crop(2, frame),frame) = smoothed;
        end
    else
        mask(1,1,frame) = NaN;
    end
end

end

