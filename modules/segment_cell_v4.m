function [ mask ] = segment_cell_v4( Imstack, varargin)
%Segment a cell from the centre of the image
% Simple Otsu thresholding
%

ImW = size(Imstack{1}{1,1},2);
ImH = size(Imstack{1}{1,1},1);
def_crop = ceil([[3/16; 13/16] * ImW; [1/4; 3/4] * ImH]);

% Keep fields and defaults up to date here:
fields = {'ellipseFitVal', 'crop', 'method', 'minSize'};
defaults = {0, def_crop, 'adaptive', 5e3};

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
mask = zeros([ size(Imstack{1}{1,1}), size(Imstack{1},1)]);

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
    threshIm = imbinarize(subIm, par.method);
    
    
    % Remove all objects smaller than minSize (in px)
    nosmallIm = bwareaopen(threshIm, par.minSize);
    
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