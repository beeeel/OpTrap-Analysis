function [ mask, fits, SegFails ] = segment_cell_v5( Imstack, varargin)
%Segment a cell from the centre of the image
% Use computer vision to segment a cell from the image given. If no values
%  for crop are given, the cell is assumed to be centrally in the image.
%
% First the image is flattened to remove noise, then sharpened with a
%  Laplacian filter to highlight edges. Active contours are used to remove
%  the background. The edges are eroded to separate noise from the circle
%  about the cell. Finally, small objects are removed and the mask is
%  returned.
%
% mask is returned as an array of size [image x, image y, num frames]. For 
%  frames where find_fails(frame) == 1, mask(1,1,frame) == NaN.
% fits has rows containing [x0, y0, a, b, theta, score]. seg_fails contains
%  0 if find_cell values for crop, centres and radius were used, 1 if
%  default values were used and 2 if the values from the previous frame
%  were used
%
% Possible input arguments:
%   k_size     - Size of flattening kernel
%   Kernel     - 'flat' or 'gaussian' kernel to denoise image
%   Lsigma     - Sigma value for Laplacian filter (edge threshold level)
%   Lalpha     - Alpha value for Laplacian filter (detail smoothing level)
%   Lbeta      - Beta value for Laplacian filter (dynamic range enhancement)
%   method     - Active contour method
%   iterations - Active contour iterations - 10 - 25 works well, depending on sc_up for find_cell_v2
%   alt_iterations - Number of iterations for when find_cell fails - should be higher than iterations
%   DiskSize   - Disk size for erosion structuring element
%   minSize    - Minimum size of objects to keep in mask (in px)
%   ellipseFitVal - How to create mask from thresholded image (0-4)
%   crop       - Image area containing the cell (opposite corners of box) - [y0; y1; x0; x1]
%   radius     - Radius from find_cell_v2
%   sc_up      - Scale factor for creating initial circular mask from find_cell's result
%   centres    - Centres vector from find_cell
%   AC_args    - Additional arguments for activecontour
%   randomize  - Randomize factor for ellipseDetection (default 2)
%
% When using a stricter crop, the number of iterations must be lower to
%prevent the contour impinging on the cell.
% AC_args: SmoothFactor defines the smoothness of the final mask - at least
%2 is good.
% ellipseFitVal: 
%0 - no fitting. 
%1 - ellipsefitfun (unknown). 
%2 - ellipseDetection and full mask of cell. 
%3 - ellipseDetection and fitted ellipse as mask.
%4 - ellipseDetection and outline of mask (as  used for fitting)

[ImH, ImW] = size(Imstack{1}{1,1});

if ImW >= 1920;     def_crop = ceil([[3/16; 13/16] * ImW; [1/4; 3/4] * ImH]);
else;               def_crop = [1; ImW; 1; ImH]; end
def_radius = 130;

% Keep fields and defaults up to date here:
fields = {'KSize', 'Lsigma', 'Lalpha', 'Lbeta', 'ACMethod', 'iterations',...
    'DiskSize', 'minSize', 'ellipseFitVal', 'crop', 'radius',...
    'alt_iterations','AC_args', 'sc_up', 'centres','randomize','Kernel',...
    'Verbose'};
defaults = {5, 250, 0.8, 2, 'Chan-Vese', 100, ...
    10,         5e3,        0,              def_crop, def_radius, ...
    200,            {},         1.2, [960; 540], 2,'gaussian',...
    false};

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
SegFails = zeros(size(Imstack{1},1),'uint8');
% Kernel for flattening
if strcmp(Par.Kernel,'flat'); Kernel = ones(Par.KSize)./(Par.KSize.^2); end
% Structuring elements for erosion
%se = strel('disk', par.DiskSize);
se2 = strel('disk',10);

for frame = 1:size(Imstack{1},1)
    % Crop out relevant area - if the input is a full stack, take the
    % relevant crop values from the crop array. If find_cell failed, crop
    % will contain NaNs, so check for those and use the default crop values
    % in that case.
    if size(Imstack{1},1) > 1 && size(Par.crop,2) > 1 && ~isnan(Par.crop(1,frame))
        Crop = Par.crop(:,frame);
        % In some cases, the crop box given will be outside the image
        % space. In these cases, replace these values with the defaults to
        % prevent indexing errors.
        Crop([1; -1; 1; -1] .* Crop < [1; -1; 1; -1] .* def_crop) = ...
            def_crop([1; -1; 1; -1] .* Crop < [1; -1; 1; -1] .* def_crop);
        Radius = Par.radius(frame);
        % Centre of circle defined relative to the field of view analysed
        Centre = Par.centres(:,frame) - Crop([1;3]);
        its = Par.iterations;
        stri = 'Using find_cell values ';
        SegFails(frame) = 0;
    elseif size(Imstack{1},1) > 1 && size(Par.crop,2) == 1
        Crop = Par.crop;
        Centre = size(Imstack{1}{frame,1})/2 - Crop([3,1])';
        Centre = Centre(end:-1:1);
        Radius = def_radius;
        its = Par.iterations;
        stri = 'Using default values ';
        SegFails(frame) = 1;
    elseif isnan(Par.crop(1,frame)) && frame > 1
        % use the same crop etc as last frame
        stri = 'Using last frame''s values ';
        SegFails(frame) = 2;
    else
        Crop = def_crop;
        Centre = [960; 540] - Crop([1;3]);
        Radius = 130;
        its = Par.alt_iterations;
        stri = 'Using default values ';
        SegFails(frame) = 1;
    end
    str = ['%s\rFrame %i:\t', stri, 'for cell location\t\t'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%imclearborder%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is where the fun begins
    subIm = Imstack{1}{frame,1}(Crop(3):Crop(4), Crop(1):Crop(2));
    % Flatten with ones kernel - maintains an amount of sharpness
    % differently to a Gaussian. I don't know if this is necessary, or
    % if it's better than a Gaussian. You're welcome to try other
    % things
    switch Par.Kernel
        case 'flat'
            % Flatten the frame using kernel above
            FiltIm = conv2(subIm, Kernel, 'same');
        case 'gaussian'
            % Flatten with gaussian filter
            FiltIm = imgaussfilt(subIm,(Par.KSize-1)/4);
    end
    
    % Sharpen image with Laplacian filter - basically a highly tunable
    % way of exaggerating the edges.
    sharpIm = locallapfilt(FiltIm, Par.Lsigma, Par.Lalpha, Par.Lbeta);
    
    % Make an initial circular mask, and use active contours to fit it to
    % the cell image.
    
    [rr, cc] = meshgrid(1:Crop(2)-Crop(1)+1,...
        1:Crop(4)-Crop(3)+1);
    preMask = (rr-Centre(1)).^2 + ...
        (cc-Centre(2)).^2 <= (Radius*Par.sc_up).^2;
    
    bw = activecontour(sharpIm-subIm, preMask, its, Par.ACMethod, Par.AC_args{:});
    
    % Erode image with disk structuring element - removes noisy bits from
    % the edge
    %erodedIm = imerode(bw,se);
    
    % Remove all objects smaller than minSize (in px)
    nosmallIm = bwareaopen(bw, Par.minSize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Par.ellipseFitVal == 1
        %Now generate Elipse paramaters and fit using Elipse fitting function
        [semiY,semiX,centreX,centreY] = ellipse_fit_fun_final(nosmallIm);
        
        mask(:,:,frame) = ellipse_mask(semiX,semiY,centreX,centreY,nosmallIm);
        
    elseif sum(Par.ellipseFitVal == [2, 3, 4]) == 1 
        %Else if ellipseFitVal ==2, use ellipseDetection method
        % Turn the mask into an outline, dilate it to reduce noise, then
        % use ellipseDetection function to find and measure the ellipse.
        % Finally, create a mask of that ellipse in the output.
        outline = imdilate(bwperim(nosmallIm),se2);
        info = regionprops(nosmallIm,'MinorAxisLength', 'MajorAxisLength','BoundingBox');
        
        % Pass a subset of the image into
        params = struct('minMajorAxis',floor(0.95 * info.MajorAxisLength),...
            'maxMajorAxis',ceil(1.05*info.MajorAxisLength),'numBest',1,...
            'randomize',Par.randomize);
        fprintf('Detecting ellipses in frame %d\n',frame)
        fits(frame,:) = ellipseDetection(outline, params);
        fits(frame,1:2) = fits(frame,1:2) + [Crop(1), Crop(3)];
        
        if Par.ellipseFitVal == 2
            mask(Crop(3):Crop(4), Crop(1):Crop(2),frame) = nosmallIm;
        elseif Par.ellipseFitVal == 3
            mask(Crop(3):Crop(4), Crop(1):Crop(2),frame)...
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
            mask(Crop(3):Crop(4), Crop(1):Crop(2), frame)...
                = outline;
        else
            error('What have you done!')
        end
    else
        mask(Crop(3):Crop(4), Crop(1):Crop(2),frame) = nosmallIm;
    end
    if Par.Verbose; fprintf([str 'segmented area = %ipx\n'],repmat(' ',1,102),frame,sum(mask(:,:,frame),'all')); end
    prog = ceil(100 * frame / size(Imstack{1},1));
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
end
fprintf('%s\r',repmat(' ',1,102))

end
