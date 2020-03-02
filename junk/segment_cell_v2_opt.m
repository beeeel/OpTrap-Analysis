function [ metric_val ] = segment_cell_v2_opt( Imstack, metric, sigma, Ethresh, Fthresh)
%Segment a cell from the centre of the image - for use with optimiser
%   Use computer vision to segment a cell from the image given. The cell is
%   assumed to be centrally in the image.
%   Input parameters are limited to image processing parameters
%   Returns a segmentation metric based on regionprops

% 4 parameters are preset - crop the middle quarter
fields = {'crop', 'minSize', 'windowSize' ...
     'ellipseFitVal', 'fails'};
defaults = {[270; 810; 480; 1440], 2000, 53, ...
    0, 0};

par = cell2struct(defaults, fields, 2);
par.sigma = sigma;
par.Ethresh = Ethresh;
par.Fthresh = Fthresh;

% Preallocate mask array - although this is not the output
mask = zeros([ size(Imstack{1}{1,1}), size(Imstack{1},1)]);
scores = zeros(size(Imstack{1},1),1);

for frame = 1:size(Imstack{1},1)
    % Crop out relevant area
    subIm = Imstack{1}{frame,1}(par.crop(1):par.crop(2), par.crop(3):par.crop(4));
    % Edge detection using canny method
    ImEdge = edge(subIm, 'canny');
    % Apply gaussian blur - the cell has lots of edges inside it
    ImFilt = imgaussfilt(double(ImEdge), par.sigma) * 10; % sigma = 6 was ok
    % Threshold to remove areas with few edges (outside cell)
    ImThresh = ImFilt > par.Ethresh;
    % Remove objects smaller than a certain size of pixels
    ImOpen = bwareaopen(ImThresh, par.minSize);
    % Apply a flattening/average filter and threshold
    kernel = ones(par.windowSize) / par.windowSize ^ 2;
    smoothed = conv2(single(ImOpen), kernel, 'same') > par.Fthresh;
    
    
    % I basically don't use anything other than ellipseFitVal = 0
    if par.ellipseFitVal == 1
        %Now generate Elipse paramaters and fit using Elipse fitting function
        [semiY,semiX,centreX,centreY] = ellipse_fit_fun_final(smoothed);
        
        mask(:,:,frame) = ellipse_mask(semiX,semiY,centreX,centreY,smoothed);
        
    elseif par.ellipseFitVal ==2
        %Else if ellipseFitVal ==2, use ellipseDetection method
        [bestFits] = ellipseDetection(ImEdge);
        
        mask(par.crop(3,frame):par.crop(4, frame), ...
            par.crop(1, frame):par.crop(2, frame),frame)...
            = ellipse_mask(majorAxis,minorAxis,centreX,centreY,bwSmoothed);
    else
        %If not fitting Ellipse, just set to previous Mask (BWSmoothed)
        mask(par.crop(1):par.crop(2), ...
            par.crop(3):par.crop(4),frame) = smoothed;
    end
    if sum(sum(mask(:,:,frame))) == 0
        scores(frame) = NaN;
    elseif strcmp(metric, 'Solidity') == 1
        props = regionprops(mask(:,:,frame), 'Solidity');
        scores(frame) = props.Solidity;
    elseif strcmp(metric,'Area') == 1
        props = regionprops(mask(:,:,frame), 'Area');
        scores(frame) = props.Area;
    else
        error('Metric not available')
    end
    
end
if sum(isnan(scores)) ~= 0
    metric_val = 10e20;
elseif strcmp(metric, 'Solidity') == 1
    error('You have not coded this metric yet');
elseif strcmp(metric,'Area') == 1
    metric_val = var(scores);
else
    error('Metric not available')
end

end