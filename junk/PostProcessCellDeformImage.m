function [info] = PostProcessCellDeform(fileName);

%4-10/aug/2017 automated post processing of .avi videos for deformation
%work Aishah Mustapha | Mina Mossayebi

%Edited: 09/May/2018 | Alex Lawrenson Made changes to the segmentation
%section to improve cell localisation

%Edited: 28/May/2018 | Alex Lawrenson Made changes to que postprocessing
%and moved all figure plotting to seperate file

% Edited: 9/May/2019 | Will Hardiman Changed input format to take ome/tiff
% stacks and included option to use alternative segmentation pipeline

% 1-read movie,
if strcmp(fileName(end-7:end), '.ome.tif')
    data = bfopen(fileName);
    Im = data{1,1}{1,1};
else
    try
        Im = imread(fileName);
    catch
        warndlg( 'File named in edit box does not appear to be a usable image file'); %in case of a wrong format of video chosen
        return
    end
end

 % Preallocate threshholdedmovie array of structures. 
 ThreshMov=...
        struct('cdata', zeros(size(Im,1), size(Im,2),1, 'uint8'),...
           'colormap', []);
        
 %preallocate an array of info structures
  info=struct('MinorAxisLength',0,'MajorAxisLength',0,'Orientation',0,'centroid',zeros(1,2), 'Perimeter',0,'Area',0,'Eccentricity',0);
  
% post process the image to be suitable to extract desired data such as
for k = 1
    
    %Read Frame

    %Temporary variable that must be set for each cell as some images are
    %at a different resolution. This should be improved by making it ratio
    %of the image size at first maybe? Or by preprocessing the images
    %quickly first to identify the min cell size?
    CellMinSize = 130;
    
    %Segment Image, can set EllipseFitVal = 1 to fit ellipses that are not
    %orientated before regionprops fits ellipses again.
    ThreshMov.cdata = segment_cell(Im, CellMinSize);
    
    %Use Regionprops to extract information and populate info structure
    
    if isequal(size(regionprops(ThreshMov.cdata,'MinorAxisLength')),[1 1])
        info.MinorAxisLength=regionprops(ThreshMov.cdata,'MinorAxisLength');
        info.MajorAxisLength=regionprops(ThreshMov.cdata,'MajorAxisLength');
        info.Orientation=regionprops(ThreshMov.cdata,'Orientation');
        info.centroid=regionprops(ThreshMov.cdata,'centroid');
        info.Perimeter=regionprops(ThreshMov.cdata,'Perimeter');
        info.Area=regionprops(ThreshMov.cdata,'Area');
        info.Eccentricity=regionprops(ThreshMov.cdata,'Eccentricity');
        
        %converts arrays of structures that we need to operate on numerically, to
        %arrays of cells so that we can use </> operators on them
        shortAxis_value=info.MinorAxisLength.MinorAxisLength;
        LongAxis_value=info.MajorAxisLength.MajorAxisLength;
        eccentricity_value=info.Eccentricity.Eccentricity;
        orientation_value=info.Orientation.Orientation;
        centroid_value_x=info.centroid.Centroid(1,1);
        centroid_value_y=info.centroid.Centroid(1,2);
        Perimeter_value=info.Perimeter.Perimeter;
        Area_value=info.Area.Area;

        %populating an array holdingthe Taylor's deformation parameter(LongestAxisLength-ShortestAxisLength)/(LongestAxisLength+ShortestAxisLength)
        TaylorParam=(LongAxis_value-shortAxis_value)/(LongAxis_value+shortAxis_value);
        info.TaylorParamater = TaylorParam;
        
        %populating 'Flatness' parameter defined as:(majorAxisLenght-MinorAxisLength)/MajorAxisLength
        flatness=(LongAxis_value-shortAxis_value)/LongAxis_value;
        info.Flatness = flatness;
        
    else
        %in case we detected multiple objects in the frame, set every parameter to 'nan' to avoid complications
        %maybe think about setting the values to the last one in order that
        %there is still continous data?
        info.MinorAxisLength=nan;
        info.MajorAxisLength=nan;
        info.Orientation=nan;
        info.centroid=nan;
        info.Perimeter=nan;
        info.Area=nan;
        info.Eccentricity=nan;
        LongAxis_value=nan;
        shortAxis_value=nan;
        orientation_value=nan;
        centroid_value_x=nan;
        centroid_value_y=nan;
        Perimeter_value=nan;
        Area_value=nan;
        
        eccentricity_value=nan;
        TaylorParam=nan;
        flatness=nan;
    end
end

end
