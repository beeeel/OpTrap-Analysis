function [info, meta] = PostProcessCellDeform_v2(Imstack, varargin)
%% info = PostProcessCellDeform_v2(Imstack, varargin)
% Automatic post processing of OME/TIFF image stacks for deformation.
% Input a cell array loaded from bfopen (or load_imstack) and inputs for
%   find_cell and/or segment_cell. Inputs in name-value pairs, e.g.:
%   PostProcessCellDeform_v2(Imstack, 'find_cell', {'Rs',[100, 200]})
% Outputs an info struct array with info about each frame. Contains NaN when
%   unable to analyse the cells.
% Contents of info struct:
%    MinorAxisLength    - of mask, measured by regionprops
%    MajorAxisLength    - ditto
%    Orientation        - ditto
%    Centroid           - ditto
%    Perimeter          - ditto
%    Area               - ditto
%    Eccentricity       - ditto
%    TaylorParameter    - of mask, calculated from Min and Maj axis lengths (from regionprops)
%    Flatness           - ditto
%    uMinorAxisLength   - of mask, measured by unwrap_cell
%    uMajorAxisLength   - ditto
%    uOrientation       - ditto
%    uTaylorParameter   - of mask, calculated from unwrap_cell
%    uOffset            - Correction to cell centre location, from unwrap_cell_v2 only
%    mask               - Binary mask, fitted by segment_cell
%    centres            - Centre of cell, found by find_cell
%    radius             - Radius of cell, ditto
%    find_fails         - Success for find_cell: 0 for success, 1 for no circles found, 2 for multiple circles found
%    seg_fails          - Initial mask source for segment_cell: 0 for find_cell, 1 for default, 2 for previous frame
%    crop               - Crop box for cell (from find_cell) [x0 x1 y0 y1]
%    filepath           - Full path for file analysed
%    NRegions           - Number of regions in mask, from regionprops
%    find_cell_v        - Version number used for find_cell
%    Seg_Cell_v         - Version number used for segment_cell
%    ellipse_fits       - Output of ellipseDetection [x0 y0 a b theta score]: centres, axis lenghts, angle and quality of fit metric
%    TotalRunTime       - Total time taken for PostProcessCellDeform to run
%
% Note: info field "fielpath" will contain the path up to the first space.
% In order for this to be the full path, directory and filenames must have
% no spaces.

% Support: william.hardiman@nottingham.ac.uk OR bill.hardiman1995@gmail.com

%4-10/aug/2017 automated post processing of .avi videos for deformation
%work Aishah Mustapha | Mina Mossayebi

%Edited: 09/May/2018 | Alex Lawrenson Made changes to the segmentation
%section to improve cell localisation

%Edited: 28/May/2018 | Alex Lawrenson Made changes to que postprocessing
%and moved all figure plotting to seperate file

%Edited: 9/May/2019 | Will Hardiman: 
% * Adapted script to accept ome/tiff image stacks (requires bioformats 
%package available from:
%https://docs.openmicroscopy.org/bio-formats/5.7.1/users/matlab/index.html)
% * Included alternative cell segmentation script
% * Changed input arguments to accept segment_cell version and parameters
%for v2
% * Reduced memory footprint by removing unused arrays (data being
%populated into info)

%Edited: 13/May/2019 | Will Hardiman:
% * Included cell locating module utilising Hough circles transform
% * Made cell_radii to iterate over cropped image stack and calculate radii
%from bright halos around cells. It didn't work.

%Edited: 16/May/2019 | Will Hardiman:
% * Forked v2 to act on cell arrays containing image data to simplify file
%loading
% * Stripped file loading and AVI support
% * Tidied struct preallocation and removed obsolete blocks to improve
%efficiency

%Edited: 17/May/2019 | Will Hardiman:
% * Restructured cell_data and masks to be put into info struct. Both are
%cleared from memory after being copied into the struct.
% * Created variables to store size of frame and number of frames, meaning
%Imstack can be cleared from memory after calling segment_cell_v2

%Edited: 28/May/2019 | Will Hardiman:
% * Created segment_cell_v3 to detect cells by a different method
% * Added switch statement to enable choosing between versions of
%segment_cell

%Edited: 04/June/2019 | Will Hardiman:
% * Added tracking for when regionprops detects multiple regions
%(NRegions field in info struct)
% * Added region selection - if multiple regions are found, take the region
%whose centroid is closest to the centre of the image (Euclidean distance)

% Edited: 06/June/2019 | Will Hardiman:
% * Added optional crop to input arguments - must have 4 rows containing 
%[x1; x2; y1; y2] where x1 and x2 are the left and right sides of the crop
%box, and y1 and y2 are the top and bottom edges. Can be either a vector,
%or a matrix with 1 column per frame.

% Edited: 07/June/2019 | Will Hardiman:
% * Created find_cell_v2 and added switch statement to choose between them.
%v2 uses a Gaussian filter, followed by a Laplacian filter to enhance the
%cell halo before using Hough circles to find the cell.

% Edited: 12/June/2019 | Will Hardiman:
% * Added optional input arguments to provide arguments to modules
%find_cell and segment_cell, including versions.

% Edited: 28/October/2019 | Will Hardiman:
% * Added segment_cell_v5 which includes a second output for the ellipse
%fitting result.

% Edited: 19/November/2019 | Will Hardiman:
% * Fixed some errors that arose when running without input arguments, and
%some that came from find_cell failing
% * Changed fails (from find_cell) to find_fails

% Edited: December/2019 | Will Hardiman:
% * Added unwrap_cell - radial sampling of the image to create an (r,theta)
%image of the cell, which can be fitted to the polar equation of an
%ellipse.
% * Added unwrap_cell_v2 - start by fitting an off-centre circle, then a
%centred ellipse. Good when find_cell isn't as accurate as one might hope.

% Edited: February/2020 | Will Hardiman:
% * Added metadata output - settings used to run each module, etc
% * Added field removal for unused fields dependent on which modules will
% be used
% * Added LineMaxima_v1 - Another way of finding the cell centre, using
% contrast enhancement to find the bright edge of the cell
% * Restructured preallocation of info and metadata to enhance readability
StartTime = tic;

% Module control parameters - these are overridden by user-defined inputs
find_cell_v = 0; % Find_cell accurately crops images after circle finding
seg_cell_v = 0; % Different versions of segment_cell use different methods to produce a mask
unwrap_cell_v = 2;
line_maxima_v = 1;

%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs
% Get file path from Imstack metadata for first frame
metadat = strsplit(Imstack{1}{1,2});
% Number of frames
N_frames = size(Imstack{1},1);
% Size of each frame
sz_frame = size(Imstack{1}{1,1});

% Process inputs - if there's an odd number greater than one.
if mod(nargin,2) == 1 && nargin > 1
    % For each even numbered input (the names)
    for argin = 1:2:nargin-1
        % If the name matches one of the expectations, put it there
        if strcmp(varargin{argin}, 'find_cell')
            disp('Setting find_cell arguments')
            find_cell_args = varargin{argin + 1};
        elseif strcmp(varargin{argin}, 'segment_cell')
            disp('Setting segment_cell arguments')
            seg_cell_args = varargin{argin + 1};
        elseif strcmp(varargin{argin}, 'unwrap_cell')
            disp('Setting unwrap_cell arguments')
            unwrap_cell_args = varargin{argin + 1};
        elseif strcmp(varargin{argin}, 'line_maxima')
            disp('Setting LineMaxima arguments')
            line_maxima_args = varargin{argin + 1};
        elseif strcmp(varargin{argin}, 'seg_cell_v')
            disp('Setting segment_cell version')
            seg_cell_v = varargin{argin + 1};
            if sum(seg_cell_v == [0, 1, 2, 3, 4, 5, 6]) ~= 1
                error(['Invalid segment_cell version ', num2str(seg_cell_v)]);
            end
        elseif strcmp(varargin{argin}, 'find_cell_v')
            disp('Setting find_cell version')
            find_cell_v = varargin{argin + 1};
            if sum(find_cell_v == [0, 1, 2]) ~= 1
                error(['Invalid find_cell version ', num2str(find_cell_v)]);
            end
        elseif strcmp(varargin{argin}, 'unwrap_cell_v')
            disp('Setting unwrap_cell version')
            unwrap_cell_v = varargin{argin + 1};
            if sum(unwrap_cell_v == [0, 1, 2]) ~= 1
                error(['Invalid unwrap_cell version ', num2str(unwrap_cell_v)]);
            end
        elseif strcmp(varargin{argin},'line_maxima_v')
            disp('Setting LineMaxima version')
            line_maxima_v = varargin{argin + 1};
            if sum(line_maxima_v == [0, 1]) ~= 1
                error(['Invalid line_maxima version ', num2str(line_maxima_v)]);
            end
        else
            % If the name isn't recognised, give an error
            error(['Unrecognised module for input argument ', num2str(argin)])
        end
    end
    % If there is an even number greater than 1, give an error
elseif nargin > 1
    error('Inputs must be in name-value pairs')
end

% Check both find_cell_args and seg_cell_args exist, and if not, populate
% them with hardcoded values
if size(whos('find_cell_args'),1) == 0 
    switch find_cell_v
        case 0
            find_cell_args = {};
        case 1
            find_cell_args = {};
        case 2
            find_cell_args = {};
    end
else
    %find_cell_v = 0;
end
if size(whos('seg_cell_args'),1) == 0
    % Segment_cell only takes input arguments for v2 onwards.
    switch seg_cell_v
        case 0
            seg_cell_args = {};
        case 2
            seg_cell_args = {'Ethresh', 0.5473, ...
            'sigma', 1, 'Fthresh', 0.5739, 'ellipseFitVal', 0};
        case 3
            seg_cell_args = {'thresh', 90, 'Lsigma', 500, 'Lalpha', 0.95,'ellipseFitVal', 0, 'threshP', 85};
            
            % Values for Lalpha and threshP have been lightly tuned from
            % defaults.
        case 4
            seg_cell_args = {};
        case 5
            seg_cell_args = {'iterations',200,'Lsigma', 0.1, 'Lalpha', 5, 'Lbeta', 10};
        case 6
            seg_cell_args = {};
    end
else
    %seg_cell_v = 0;
end
if size(whos('unwrap_cell_args'),1) == 0 
    switch unwrap_cell_v
        case 0 
            unwrap_cell_args = {};
        case 1
            unwrap_cell_args = {};
        case 2
            unwrap_cell_args = {};
    end
else
    %unwrap_cell_v = 0;
end
if isempty(whos('line_maxima_args'))
    switch line_maxima_v
        case 0 
            line_maxima_args = {};
        case 1
            line_maxima_args = {};
    end
else
    %line_maxima_v = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preallocate an array of info structures

% The order of fields in {regionFields, fields} must be the same as the
% order of defaults in values. Hopefully the linebreaks help?
RegionFields = {'MinorAxisLength', 'MajorAxisLength', 'Orientation', ...
    'Centroid', 'Perimeter', 'Area', 'Eccentricity'};
RegionValues = {0,        0,      0, ...
    zeros(1,2),         0,      0,      0,};

% Concatenating cell arrays with [] is quicker than deconstructing with {:}
% and reconstructing inside {...}
SegFields = [RegionFields, 'mask', 'seg_fails','ellipse_fits', 'TaylorParameter', 'Flatness'];
SegValues = [RegionValues, false(sz_frame), uint8(0), zeros(6,1), 0, 0];

FindFields = {'centres', 'radius', 'find_fails', 'crop'};
FindValues = {[0;0], 0, uint8(0), [1;sz_frame(2);1;sz_frame(1)]};

UnwrapFields = {'uMajorAxisLength','uMinorAxisLength',...
    'uOrientation','uTaylorParameter','uFlatness','uOffset','uFitErrs'};
UnwrapValues = {0, 0, ...
    0, 0, 0, zeros(3,1), zeros(3,1)};
InfoFields = [FindFields, UnwrapFields, SegFields, 'filepath', 'mCentres'];
InfoValues = [FindValues, UnwrapValues, SegValues, metadat{1}, [0;0]];

% I'm leaving this here in case I broke something
% InfoFields = [SegFields, 'TaylorParameter', 'Flatness', 'mask', ...
%     'centres', 'radius', 'find_fails', 'seg_fails', ...
%     'crop', 'filepath', 'uMajorAxisLength','uMinorAxisLength',...
%     'uOrientation','uTaylorParameter','uFlatness','uOffset','ellipse_fits',...
%     'mCentres'];
% InfoValues = [SegValues, 0,  0,  false(sz_frame),    ...
%     [0; 0], 0,  uint8(0), uint8(0), ...
%     [1; sz_frame(2); 1; sz_frame(1)], metadat{1}, 0, 0, ...
%     0, 0, 0, zeros(3,1), zeros(6,1),...
%     [0; 0]];

MetaFields = {'filepath', 'N_Frames','Frame_size', 'TotalRunTime', ...
    'Find_cell_time', 'Unwrap_cell_time', 'Segment_cell_time', 'Line_maxima_time',...
    'Find_cell_args', 'find_cell_v', 'Segment_cell_args', 'seg_cell_v', ...
    'Unwrap_cell_args', 'unwrap_cell_v', 'Line_maxima_args', 'line_maxima_v'};
MetaValues = {metadat{1}, N_frames, sz_frame,0, ...
    0, 0, 0, 0,...
    find_cell_args, find_cell_v, seg_cell_args, seg_cell_v, ...
    unwrap_cell_args, unwrap_cell_v, line_maxima_args, line_maxima_v};

% Create one struct for info and one for metadata
info(1:N_frames) = cell2struct(InfoValues, InfoFields, 2);
meta = cell2struct(MetaValues, MetaFields, 2);

% Tidy it to remove excess fields
if seg_cell_v == 0
    info = rmfield(info,SegFields );
    meta = rmfield(meta, {'Segment_cell_time','Segment_cell_args'});
end
if find_cell_v == 0
    info = rmfield(info, FindFields);
    meta = rmfield(meta, {'Find_cell_time','Find_cell_args'});
end
if unwrap_cell_v == 0
    info = rmfield(info, UnwrapFields);
    meta = rmfield(meta, {'Unwrap_cell_time','Unwrap_cell_args'});
end
if line_maxima_v == 0
    info = rmfield(info, 'mCentres');
    meta = rmfield(meta, {'Line_maxima_time','Line_maxima_args'});
end

%% Run find_cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if find_cell_v ~= 0
    before = toc(StartTime);
    % 1-Find cells,
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Finding cells');
    disp(['Using find_cell version ',num2str(find_cell_v)]);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    % Starting cell finding - returns a struct array with info on the best
    % circle in each image, tracking of when it failed, and a crop array
    % showing where to crop for the best circle in each image.
    switch find_cell_v
        case 1
            cell_dat = find_cell(Imstack, find_cell_args{:});
        case 2
            cell_dat = find_cell_v2(Imstack, find_cell_args{:});
    end
    
    % Transfer the info from cell_dat struct into info struct
    InfoFields = fieldnames(cell_dat);
    for frame = 1:N_frames
        for f_no = 1:numel(InfoFields)
            info(frame).(InfoFields{f_no}) = cell_dat(frame).(InfoFields{f_no});
        end
    end
    clear cell_dat
    meta.Find_cell_time = toc(StartTime) - before;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Found cells');
    fprintf('%g s elapsed\n',meta.Find_cell_time)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
else
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Skipping find_cell');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%% Use LineMaxima method
if line_maxima_v ~= 0
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Finding cells with line maxima');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    before = toc(StartTime);
    switch line_maxima_v
        case 1
            Centres = LineMaxima_v1(Imstack,line_maxima_args{:});
        otherwise
            error('huhnknown line_maxima_v')
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Found cells with line maxima')
    fprintf('%g s elapsed\n',toc(StartTime))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Parsing output')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    for frame = 1:N_frames
        info(frame).mCentres = Centres(:,frame);
    end
    meta.Line_maxima_time = toc(StartTime) - before;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Found cells and finished')
    fprintf('%g s elapsed\n',toc(StartTime))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
end

%% Unwrap cells
if unwrap_cell_v ~= 0
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Unwrapping cells');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    before = toc(StartTime);
    if line_maxima_v == 0
        Centres = [info.centres];
        Radii = [info.radius];
    else
        Centres = [info.mCentres];
        Radii = repmat(min(sz_frame)/2,1,N_frames);
    end
    switch unwrap_cell_v
        case 1
            UnwrapFits = unwrap_cell_v1(Imstack, Centres, Radii, unwrap_cell_args{:});
        case 2
            [UnwrapFits, ~, ~, ~, UnwrapOffset, FitErrs] = unwrap_cell_v2(Imstack, Centres, Radii, unwrap_cell_args{:});
        otherwise
            error('huh')
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Unwrapped cells')
    fprintf('%g s elapsed\n',toc(StartTime))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Parsing output')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    for frame = 1:N_frames
        info(frame).uMajorAxisLength = UnwrapFits(1,frame);
        info(frame).uMinorAxisLength = UnwrapFits(2,frame);
        info(frame).uOrientation = UnwrapFits(3,frame);
        info(frame).uOffset = UnwrapOffset(:,frame);
        info(frame).uFitErrs = FitErrs(1:3,frame);
        
        % Taylor's deformation parameter : (LongestAxisLength-ShortestAxisLength)/(LongestAxisLength+ShortestAxisLength)
        info(frame).uTaylorParameter = ( UnwrapFits(1,frame) - ...
            UnwrapFits(2,frame)) / ( UnwrapFits(1,frame) ...
            + UnwrapFits(2,frame));
        
        % Flatness : (majorAxisLength-MinorAxisLength)/MajorAxisLength
        info(frame).uFlatness = ( UnwrapFits(1,frame) - ...
            UnwrapFits(2,frame)) /  UnwrapFits(1,frame);
    end
    meta.Unwrap_cell_time = toc(StartTime) - before;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Unwrapped cells and finished')
    fprintf('%g s elapsed\n',toc(StartTime))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    clear Centres Radii
end

if seg_cell_v ~= 0
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Masking cells');
    fprintf('Using segment_cell_v%d\n', seg_cell_v);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    before = toc(StartTime);
    %% Run segment_cell
    % Segment cell to get masks array - segment_cell_v2 and newer iterate over
    % image stack.
    switch seg_cell_v
        case 1
            % V1: Use edge detection, dilation and erosion
            % I've not had success using this
            masks = zeros([sz_frame, N_frames],'uint8');
            for frame = 1:N_frames
                masks(:,:,frame) = segment_cell(Imstack{1}{frame,1});
            end
        case 2
            % V2: Use Canny edge detection, followed by thresholding and filtering
            % This only works for very nice looking cells
            masks = segment_cell_v2(Imstack, 'crop', [info.crop],  ...
                'fails', [info.find_fails], seg_cell_args{:});
        case 3
            % V3: Use Laplacian filtering for edge enhancement, thresholding and
            % dilation.
            % This is more robust than v2.
            masks = segment_cell_v3(Imstack, 'crop', [info.crop], seg_cell_args{:});
        case 4
            masks = segment_cell_v4(Imstack);
        case 5
            if find_cell_v ~=0
                [masks, fits, SegFails] = segment_cell_v5(Imstack, 'crop', [info.crop], ...
                    'radius', [info.radius], 'fails', [info.find_fails],...
                    'centres', [info.centres], seg_cell_args{:});
            else
                [masks, fits, SegFails] = segment_cell_v5(Imstack, seg_cell_args{:});
            end
        case 6
            [masks, fits, SegFails] = segment_cell_v6(Imstack, seg_cell_args{:});
        otherwise
            error('Unknown segment_cell version %d', seg_cell_v);
    end
    
    
    %% Put masks into info struct
    for frame = 1:N_frames
        info(frame).mask = masks(:,:,frame);
        info(frame).seg_fails = SegFails(frame);
        if seg_cell_v == 5
            info(frame).ellipse_fits = fits(frame, :)';
        end
    end
    clear masks
    clear Imstack
    clear fits
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Masked cells');
    fprintf('%g s elapsed\n',toc(StartTime))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Measuring cells');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    %% 4- Get region properties
    for frame = 1 : N_frames
        % (optional:) Only look at the cells that were found - if you change
        % the number below to 1, it will skip frames where find_cell failed.
        if 1% info(frame).find_fails ~= -1
            % Use Regionprops to extract information and populate info structure
            props = regionprops(imfill(info(frame).mask,'holes'),RegionFields);
            % Put the properties into the info struct - iterate over the list
            % of field names.
            % If the mask is empty, props is a 0x1 struct array, which creates
            % errors for dot indexing. In this case, put NaNs into info.
            
            % If regionprops finds multiple regions in the image (e.g.: other
            % cells in FoV), the distance to the centre for each centroid is
            % calculated. The object closest to the centre is chosen.
            Names = fieldnames(props);
            if size(props,1) > 0
                % Store the number of found objects
                info(frame).NRegions = size(props,1);
                % If only one object is found, just take it
                if size(props,1) == 1
                    for field = 1:numel(Names)
                        info(frame).(Names{field}) = props.(Names{field});
                    end
                else
                    % If multiple are found, calculate the distances
                    dists = zeros(size(props));
                    for idx = 1:size(props,1)
                        dists(idx) = sum((props(idx).Centroid - [960, 540]).^2);
                    end
                    % Choose the closest to the centre and take it
                    [~, idx] = min(dists);
                    for field = 1:numel(Names)
                        info(frame).(Names{field}) = props(idx).(Names{field});
                    end
                end
            else
                % If none are found, store NaNs
                for field = 1:numel(Names)
                    if ~strcmp(Names{field},'Centroid')
                        info(frame).(Names{field}) = NaN;
                    else
                        info(frame).(Names{field}) = [NaN, NaN];
                    end
                end
            end
            
            % Taylor's deformation parameter : (LongestAxisLength-ShortestAxisLength)/(LongestAxisLength+ShortestAxisLength)
            info(frame).TaylorParameter = (info(frame).MajorAxisLength - ...
                info(frame).MinorAxisLength) / (info(frame).MajorAxisLength ...
                + info(frame).MinorAxisLength);
            
            % Flatness : (majorAxisLength-MinorAxisLength)/MajorAxisLength
            info(frame).Flatness = (info(frame).MajorAxisLength - ...
                info(frame).MinorAxisLength) / (info(frame).MajorAxisLength);
            
        else
            % When no circle is detected, there's no mask, and so nothing we
            % can do :(
            % This isn't strictly true - segment_cell_v3 is good enough to not
            % need the crop.
            info(frame).MinorAxisLength=nan;
            info(frame).MajorAxisLength=nan;
            info(frame).Orientation=nan;
            info(frame).centroid=nan;
            info(frame).Perimeter=nan;
            info(frame).Area=nan;
            info(frame).Eccentricity=nan;
            info(frame).TaylorParameter=nan;
            info(frame).Flatness=nan;
        end
        prog = ceil(100 * frame / N_frames);
        fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])

    end
    fprintf('%s\r',repmat(' ',1,104))
    meta.Segment_cell_time = toc(StartTime) - before;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Measured cells and finished');
    fprintf('%g s elapsed\n',toc(StartTime))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
else
    disp('Skipping segment_cell')
end
meta.TotalRunTime = toc(StartTime) - StartTime;
end
