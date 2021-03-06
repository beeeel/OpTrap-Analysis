function [ mask, fits, SegFails ] = segment_cell_v6( Imstack, varargin)
%Segment a cell from the centre of the image
% Use computer vision to segment a cell from the image given. Does not care
% whereabouts the cell is in the image
%
% method == 'imclearborder'
% Use imclearborder to find light elements (mostly the cell), then remove
% small objects and fill holes to arrive at a segmented cell.
%
% method == 'lazysnapping' (default)
% Use superpixels and lazysnapping to segment the cell, assuming the
% brightest and darkest 'Percent'% of pixels are within the cell
%
% mask is returned as an array of size [image x, image y, num frames].
% fits has rows containing [x0, y0, a, b, theta, score]. SegFails contains
% 1 if the mask is empty, 0 otherwise.
%
% Possible input arguments:
%   minSize     - Minimum size of objects to keep in mask (in px)
%   ellipseFitVal - How to create mask from thresholded image (0-4)
%   randomize   - Randomize factor for ellipseDetection (default 2)
%   method      - Which CV approach to use
%   NSuperPixels- Number of superpixels to create (N/A to imclearborder)
%   Percent     - Bottom % and top % of pixels assumed to be within the cell
%   Kernel      - Preprocessing kernel used
%   KSize       - Size of kernel used for preprocessing
%
% ellipseFitVal: 
%0 - no fitting. 
%1 - ellipsefitfun (unknown). 
%2 - ellipseDetection and full mask of cell. 
%3 - ellipseDetection and fitted ellipse as mask.
%4 - ellipseDetection and outline of mask (as  used for fitting)

% Keep fields and defaults up to date here:
fields = {'minSize', 'ellipseFitVal', 'randomize','method', 'NSuperPixels',...
    'Percent', 'Kernel','KSize'};
defaults = {1000, 0, 2,'lazysnapping', 512,...
    2, 'gaussian', 5};
            
Par = cell2struct(defaults, fields,2);

% Parse inputs and create struct with parameters given
if ~isempty(varargin)
    if mod(nargin,2) == 0 % Imstack counts towards nargin
        error('Please supply arguments in name-value pairs');
    end
    % For each field, display it depending on its size and contents
    disp('Input arguments for segment_cell_v5:')
    for field = 1:(nargin - 1)/2
        if ischar(varargin{2*field})
            disp([varargin{2*field -1} ' = ' varargin{2*field}])
        elseif size(varargin{2 * field}, 2) < 4 && ...
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

% Preallocate output arrays
mask = false([ size(Imstack{1}{1,1}), size(Imstack{1},1)]);
fits = zeros(size(Imstack{1},1), 6);
SegFails = zeros(size(Imstack{1},1), 1);

% Structuring elements for erosion
se2 = strel('disk',10);
% Flattening kernel
if strcmp(Par.Kernel,'flat'); Kernel = ones(Par.KSize)./(Par.KSize.^2); end
% Size of images
[ImH, ImW] = size(Imstack{1}{1,1});

for frame = 1:size(Imstack{1},1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is where the fun begins
    
    switch Par.method
        case 'imclearborder'
            % Remove dark things connected to image border
            Rough = imclearborder(Imstack{1}{frame,1})~=0;

            % Remove small objects
            NoSmall = bwareaopen(Rough, Par.minSize);

            % Fill holes
            NoHoles = imfill(NoSmall, 'holes');

        case 'lazysnapping'    
            switch Par.Kernel
                case 'flat'
                    % Flatten the frame using kernel above
                    FiltIm = conv2(Imstack{1}{1,1}, Kernel, 'same');
                case 'gaussian'
                    % Flatten with gaussian filter
                    FiltIm = imgaussfilt(Imstack{1}{1,1},(Par.KSize-1)/4);
            end
            % Create superpixels
            Labels = superpixels(FiltIm, Par.NSuperPixels);
            
            % Foreground mask from the brightest and darkest bits of the
            % image (edge and middle of cell, normally)
            ForeMask = (Imstack{1}{frame,1} < prctile(Imstack{1}{frame,1},Par.Percent,'all')) + ...
                (Imstack{1}{frame,1} > prctile(Imstack{1}{frame,1},100-Par.Percent,'all'));
            
            % Background mask from the edge pixels of the image
            BackMask = zeros(ImH, ImW);
            idx = sub2ind([ImH, ImW],[ones(1,ImW-1),1:ImH-1, ImH * ones(1,ImW-1), ImH:-1:2],...
                [1:ImW-1, ImW * ones(1,ImH-1), ImW:-1:2, ones(1,ImH-1)]);
            BackMask(idx) = 1;
            
            % Lazysnapping to do the segmentation
            BW = lazysnapping(FiltIm, Labels, ForeMask, BackMask);
            
            % Remove small objects
            NoSmall = bwareaopen(BW, Par.minSize);
            
            % Remove stuff touching the edge
            NoHoles = imclearborder(NoSmall);
    end
    if sum(NoHoles,'all') == 0
        SegFails(frame) = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Par.ellipseFitVal == 1
        %Now generate Elipse paramaters and fit using Elipse fitting function
        [semiY,semiX,centreX,centreY] = ellipse_fit_fun_final(NoHoles);
        
        mask(:,:,frame) = ellipse_mask(semiX,semiY,centreX,centreY,NoHoles);
        
    elseif sum(Par.ellipseFitVal == [2, 3, 4]) == 1 
        %Else if ellipseFitVal ==2, use ellipseDetection method
        % Turn the mask into an outline, dilate it to reduce noise, then
        % use ellipseDetection function to find and measure the ellipse.
        % Finally, create a mask of that ellipse in the output.
        outline = imdilate(bwperim(NoHoles),se2);
        info = regionprops(NoHoles,'MinorAxisLength', 'MajorAxisLength');
        
        % Pass a subset of the image into
        params = struct('minMajorAxis',floor(0.95 * info.MajorAxisLength),...
            'maxMajorAxis',ceil(1.05*info.MajorAxisLength),'numBest',1,...
            'randomize',Par.randomize);
        fprintf('Detecting ellipses in frame %d\n',frame)
        fits(frame,:) = ellipseDetection(outline, params);
        fits(frame,1:2) = fits(frame,1:2) + [crop(1), crop(3)];
        
        if Par.ellipseFitVal == 2
            mask(crop(3):crop(4), crop(1):crop(2),frame) = NoHoles;
        elseif Par.ellipseFitVal == 3
            mask(crop(3):crop(4), crop(1):crop(2),frame)...
                ... % Below defines the mask for an ellipse ((x-x0)/a).^2 + ((y-y0)/b).^2 <= 1
                = sqrt(((rr-fits(frame,1))./fits(frame,3)).^2 + ...
                ((cc-fits(frame,2))./fits(frame,4)).^2) <= 1;
            % The ellipse above is always oriented with major axis along x, so
            % it needs to be rotated to match the cell orientation. There's
            % a better way to do it with a coordinate system change and
            % some trignometry, but I rarely use this setting
            mask(:,:,frame) = imrotate(mask(:,:,frame),fits(frame,5),'nearest','crop');
        elseif Par.ellipseFitVal == 4
            % Set the mask to the outline used for ellise fitting - this
            % can be used as a debug mode to check what the ellipse is
            % being fitted to.
            mask(crop(3):crop(4), crop(1):crop(2), frame)...
                = outline;
        else
            error('What have you done!')
        end
    else
        mask(:,:,frame) = NoHoles;
    end
    fprintf('%s\rFrame %i \tsegmented area = %ipx\n',repmat(' ',1,104),frame,sum(mask(:,:,frame),'all'))
    prog = ceil(100 * frame / size(Imstack{1},1));
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
end
fprintf('%s\r',repmat(' ',1,104))

end
