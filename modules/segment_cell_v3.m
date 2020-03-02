function [ mask ] = segment_cell_v3( Imstack, varargin)
%Segment a cell from the centre of the image
% Use computer vision to segment a cell from the image given. If no values
%  for crop are given, the cell is assumed to be centrally in the image.
% First the image is flattened to remove noise, then sharpened with a
%  Laplacian filter to highlight edges, and thresholded to remove
%  background. The edges are dilated and eroded to try and create a circle
%  about the cell. Finally, small objects are removed and the mask is
%  returned. 
% Mask is returned as an array of size [image x, image y, num frames].
% For frames where fails(frame) == 1, mask(1,1,frame) == NaN.
% Possible input arguments:
%   k_size     - Size of flattening kernel
%   Lsigma     - Sigma value for Laplacian filter (edge threshold level)
%   Lalpha     - Alpha value for Laplacian filter (detail smoothing level)
%   Lbeta      - Beta value for Laplacian filter (dynamic range enhancement)
%   thresh     - Percentile to threshold filtered image
%   LineLength - Line length for dilation structuring elements
%   diamondsize- Size of diamond for erosion structuring element
%   minSize    - Minimum size of objects to keep in mask (in px)
%   ellipseFitVal - How to create mask from thresholded image (legacy)
%   crop       - Image area containing the cell (opposite corners of box)
%
% When using a stricter crop, the percentile threshold, thresh, will need to
% be lower as there is less chance of other cells in the FoV.

ImW = size(Imstack{1}{1,1},2);
ImH = size(Imstack{1}{1,1},1);
def_crop = ceil([[3/16; 13/16] * ImW; [1/4; 3/4] * ImH]);

% Keep fields and defaults up to date here:
fields = {'k_size', 'Lsigma', 'Lalpha', 'Lbeta', 'method', ...
    'LineLength', 'diamondSize', 'minSize', 'ellipseFitVal', 'crop'};
defaults = {5, 250, 0.8, 2, 'global', ...
    4, 1, 5e3, 0, def_crop};

par = cell2struct(defaults, fields,2);

% Parse inputs and create struct with parameters given
if nargin > 1
    if mod(nargin,2) == 0 % Imstack counts towards nargin
        error('Please supply arguments in name-value pairs');
    end
    % For each field, display it depending on its size
    disp('Input arguments for segment_cell_v3:')
    for field = 1:(nargin - 1)/2
        if size(varargin{2 * field}, 2) < 4 && ...
                size(varargin{2 * field},1) == 1
            disp([varargin{2 * field - 1}, ' = ', num2str(varargin{2 * field})]);
        else
            disp([varargin{2 * field - 1}, ' = size[', ...
                size(varargin{2 * field},1), size(varargin{2 * field},2)...
                , ']']);
        end
        % And put it into the par struct
        par.(varargin{2*field - 1}) = varargin{2*field};
    end
end

% For each name in fields, ensure par has that field
for f = 1:size(fields,2)
    if ~isfield(par, fields{f})
        par.(fields{f}) = defaults{f};
    end
end

% Preallocate output array
mask = zeros([ size(Imstack{1}{1,1}), size(Imstack{1},1)],'uint8');
% Kernel for flattening
kernel = ones(par.k_size)./par.k_size^2;
% Structuring elements for dilation
se90 = strel('line', par.LineLength, 90);
se0 = strel('line', par.LineLength, 0);
% Structuring element for erosion
seD = strel('diamond',par.diamondSize);

for frame = 1:size(Imstack{1},1)
    % Crop out relevant area - if the input is a full stack, take the
    % relevant crop values from the crop array. If find_cell failed, crop
    % will contain NaNs, so check for those and use the default crop values
    % in that case. 
    if size(Imstack{1},1) > 1 && size(par.crop,2) > 1 && ~isnan(par.crop(1,frame))
        subIm = Imstack{1}{frame,1}(par.crop(3,frame):par.crop(4, frame), par.crop(1, frame):par.crop(2, frame));
    else
        subIm = Imstack{1}{frame,1}(def_crop(3):def_crop(4), def_crop(1):def_crop(2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is where the fun begins
    
    % Flatten with ones kernel - maintains an amount of sharpness
    % differently to a Gaussian. I don't know if this is necessary, or
    % if it's better than a Gaussian. You're welcome to try other
    % things
    flatIm = uint16(conv2(subIm, kernel, 'same'));
    
    % Sharpen image with Laplacian filter - basically a highly tunable
    % way of finding the edges.
    sharpIm = locallapfilt(flatIm, par.Lsigma, par.Lalpha, par.Lbeta);
    
    threshIm = imbinarize(sharpIm, par.method);
    % Dilate image with line structuring elements
    dilatedIm = imdilate(threshIm, [se90 se0]);
    dilatedIm = imdilate(dilatedIm, [se90 se0]);
    dilatedIm = imdilate(dilatedIm, [se90 se0]);
    
    % Erode image with diamond structuring element
    erodedIm = imerode(dilatedIm,seD);
    erodedIm = imerode(erodedIm,seD);
    erodedIm = imerode(erodedIm,seD);
    erodedIm = imerode(erodedIm,seD);
    
    % Remove all objects smaller than minSize (in px)
    nosmallIm = bwareaopen(erodedIm, par.minSize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if par.ellipseFitVal == 1
        %Now generate Elipse paramaters and fit using Elipse fitting function
        [semiY,semiX,centreX,centreY] = ellipse_fit_fun_final(nosmallIm);
        
        mask(:,:,frame) = ellipse_mask(semiX,semiY,centreX,centreY,nosmallIm);
        
    elseif par.ellipseFitVal ==2
        %Else if ellipseFitVal ==2, use ellipseDetection method
        % I don't think this actually works
        [bestFits] = ellipseDetection(nosmallIm);
        
        mask(par.crop(3,frame):par.crop(4, frame), ...
            par.crop(1, frame):par.crop(2, frame),frame)...
            = ellipse_mask(majorAxis,minorAxis,centreX,centreY,bwSmoothed);
    else
        if size(par.crop,2) == size(Imstack{1},1) && ~isnan(par.crop(1,frame))
            %If not fitting Ellipse, just set to previous mask
            mask(par.crop(3,frame):par.crop(4, frame), ...
                par.crop(1, frame):par.crop(2, frame),frame) = nosmallIm;
        else
            mask(def_crop(3):def_crop(4), def_crop(1):def_crop(2), ...
                frame) = nosmallIm;
        end
    end
end

end
