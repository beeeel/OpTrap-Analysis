function [info] = PostProcessCellDeform(fileName, varargin)
%Automatic post processing of .avi videos or OME/TIFF image stacks for
%deformation
% Second input parameter refers to which version of segment_cell to use -
% v1 is (Alex?)'s code, v2 is Will's version

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
tic
% 0 - Input argument parsing
for arg = 1:nargin-1
    if isstruct(varargin{arg})
        params = varargin{arg};
    elseif isnumeric(varargin{arg}) && sum(size(varargin{arg}) == [1,1]) == 2
        v_num = varargin{arg};
    else
        warndlg('Unrecognised input argument - should be parameters struct or segment_cell version (number)');
    end
end

% 1-read movie,

% If the file is an OME/TIFF stack, use bfopen
if strcmp(fileName(end-7:end), '.ome.tif')
    try
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Loading image stack');
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        % OME/TIFF data is loaded into a struct. The actual image data is
        % within data{1}{frame_num, 1}
        data = bfopen(fileName);
        
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Loaded image stack');
        toc
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    catch
        warndlg( 'File does not appear to be a usable OME/TIFF file');
        return
    end
    % Create a video struct so other parts don't break
    video = struct('NumberOfFrames',size(data{1},1), ...
        'Height', size(data{1}{1,1},1), 'Width', size(data{1}{1,1},2));
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Finding cells');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    % Starting cell finding
    cell_dat = find_cell(data);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Found cells');
    toc
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
        
% Otherwise it should be an AVI - use the VideoReader
else
    try
        video = VideoReader(fileName);
    catch
        warndlg( 'File named in edit box does not appear to be a usable AVI file'); %in case of a wrong format of video chosen
        return
    end
end


% Preallocate threshholdedmovie array of structures.
ThreshMov(1:video.NumberOfFrames) = ...
    struct('cdata', zeros(video.Height, video.Width,1, 'uint8'),...
    'colormap', []);

%preallocate an array of info structures
info(1:video.NumberOfFrames) = ...
    struct('MinorAxisLength', 0, 'MajorAxisLength', 0, ...
    'Orientation', 0, 'centroid', zeros(1,2), 'Perimeter', 0, ...
    'Area', 0, 'Eccentricity', 0);

% post process the image to be suitable to extr700, 1200, 100, 600act desired data such as
for k = 1 : video.NumberOfFrames
    
    %Read Frame
    if strcmp(fileName(end-7:end), '.ome.tif')
        Im = data{1}{k, 1};
    else
        Im = rgb2gray(read(video, k));
    end
    
    %Temporary variable that must be set for each cell as some images are
    %at a different resolution. This should be improved by making it ratio
    %of the image size at first maybe? Or by preprocessing the images
    %quickly first to identify the min cell size?
    CellMinSize = 13000;
    
    %Segment Image, can set EllipseFitVal = 1 to fit ellipses that are not
    %orientated before regionprops fits ellipses again.
    if v_num == 1
        ThreshMov(k).cdata = segment_cell(Im, CellMinSize);
    elseif v_num == 2
        ThreshMov(k).cdata = segment_cell_v2(Im, params);
    end
    

    %Use Regionprops to extract information and populate info structure
    
    if isequal(size(regionprops(ThreshMov(k).cdata,'MinorAxisLength')),[1 1])
        info(k).MinorAxisLength=regionprops(ThreshMov(k).cdata,'MinorAxisLength');
        info(k).MajorAxisLength=regionprops(ThreshMov(k).cdata,'MajorAxisLength');
        info(k).Orientation=regionprops(ThreshMov(k).cdata,'Orientation');
        info(k).centroid=regionprops(ThreshMov(k).cdata,'centroid');
        info(k).Perimeter=regionprops(ThreshMov(k).cdata,'Perimeter');
        info(k).Area=regionprops(ThreshMov(k).cdata,'Area');
        info(k).Eccentricity=regionprops(ThreshMov(k).cdata,'Eccentricity');
        
        %converts arrays of structures that we need to operate on numerically, to
        %arrays of cells so that we can use </> operators on them
        shortAxis_value=info(k).MinorAxisLength.MinorAxisLength;
        LongAxis_value=info(k).MajorAxisLength.MajorAxisLength;

        %populating an array holdingthe Taylor's deformation parameter(LongestAxisLength-ShortestAxisLength)/(LongestAxisLength+ShortestAxisLength)
        TaylorParam=(LongAxis_value-shortAxis_value)/(LongAxis_value+shortAxis_value);
        info(k).TaylorParamater = TaylorParam;
        
        %populating 'Flatness' parameter defined as:(majorAxisLenght-MinorAxisLength)/MajorAxisLength
        flatness=(LongAxis_value-shortAxis_value)/LongAxis_value;
        info(k).Flatness = flatness;
        
    else
        %in case we detected multiple objects in the frame, set every parameter to 'nan' to avoid complications
        %maybe think about setting the values to the last one in order that
        %there is still continous data?
        info(k).MinorAxisLength=nan;
        info(k).MajorAxisLength=nan;
        info(k).Orientation=nan;
        info(k).centroid=nan;
        info(k).Perimeter=nan;
        info(k).Area=nan;
        info(k).Eccentricity=nan;
        info(k).TaylorParameter=nan;
        info(k).Flatness=nan;
    end
end


end