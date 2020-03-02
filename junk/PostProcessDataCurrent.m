%% Script for plotting data already collected from postprocessing script
% Reads saved data from processed Cell trapping videos (in .mat form) and
% then allows plotting of the different fields. Enabling Code folding under 
% preferences is recommended. Developed for use with PostProcessCellDeform 
% written by the same authors. mat file selected must contain a struct with
% several fields
% 
% Last Updated: 01/Jun/2018
% Authors: Aishah Mustapha | Mina Mossayebi | Alex Lawrenson
% Email: alexlawrenson2@gmail.com
%
% TO DO:
%   - consider making this into a function that can operate on several data
%   files
%   - could get it to do all of the statistics and return a figure showing
%   the final variables required. Taylor Parameter etc.
%   - Though would need to make it take account of the different file types
%   ie. some were for single movements, others movements back and forth,
%   others several repeat movements.
%   - Need to possible integrate with Brownian motion file and then
%   everything can be done in one run
%   - currently relies on field names predefined when running
%   PostProcessCellDeform. Need to change this to make it dynamic maybe?

%% Reading and processing saved data files
clear all; close all; clc;

%User selects file to load data from
[chosenFile, chosenpath] = uigetfile('*.mat', 'Select data to process');
fileName = fullfile(chosenpath, chosenFile);

%load struct and asign data to CellDeformInfo
tempStruct = load(fileName);
mYfieldName = fieldnames(tempStruct); %have to find the field name first
CellDeformInfo = tempStruct.(mYfieldName{1}); %then access this field with new variable

nFrames = length(CellDeformInfo);

%Now read data for every frame
for k = 1 : nFrames

    %if statement to insure only finite values are read, needed for nans
    if(isstruct([CellDeformInfo(k).MinorAxisLength]))
        %converts arrays of structures that we need to operate on numerically, to
        %arrays of cells so that we can use </> operators on them
        LongAxis_value(k)=CellDeformInfo(k).MajorAxisLength.MajorAxisLength;
        shortAxis_value(k)=CellDeformInfo(k).MinorAxisLength.MinorAxisLength;
        
        eccentricity_value(k)=CellDeformInfo(k).Eccentricity.Eccentricity;
        orientation_value(k)=CellDeformInfo(k).Orientation.Orientation;
        Perimeter_value(k)=CellDeformInfo(k).Perimeter.Perimeter;
        Area_value(k)=CellDeformInfo(k).Area.Area;
        
        %populating an array holdingthe Taylor's deformation parameter(LongestAxisLength-ShortestAxisLength)/(LongestAxisLength+ShortestAxisLength)
        TaylorParam(k)=(LongAxis_value(k)-shortAxis_value(k))/(LongAxis_value(k)+shortAxis_value(k));
        
        %populating 'Flatness' parameter defined as:(majorAxisLenght-MinorAxisLength)/MajorAxisLength
        flatness(k)=(LongAxis_value(k)-shortAxis_value(k))/LongAxis_value(k);
        
    else
        %if struct contains nans then set these seperately
        LongAxis_value(k)=nan;
        shortAxis_value(k)=nan;
        orientation_value(k)=nan;
        centroid_value_x(k)=nan;
        centroid_value_y(k)=nan;
        Perimeter_value(k)=nan;
        Area_value(k)=nan;
        eccentricity_value(k)=nan;
        TaylorParam(k)=nan;
        flatness(k)=nan;
    end
end

    %finding the maximum and minimum of major and minor axis lengths.
    [maxMajorAxisLength,I_maxMajorAxisLength]=max(LongAxis_value);
    [minMajorAxisLength,I_minMajorAxisLength]=min(LongAxis_value);
    [maxMinorAxisLength,I_maxMinorAxisLength]=max(shortAxis_value);
    [minMinorAxisLength,I_minMinorAxisLength]=min(shortAxis_value);
    
    %Taylor Param
    [maxTaylor,I_maxTaylor]= max(TaylorParam);
    [minTaylor,I_minTaylor]= min(TaylorParam);
    
    %Eccentricity
    [maxEccentricity,I_MaxEccentricity]=max(eccentricity_value);
    [minEccentricity,I_MinEccentricity]=min(eccentricity_value);
    
    %flattness
    [maxFlatness,I_maxFlatness]= max(TaylorParam);
    [minFlatness,I_minFlatness]= min(flatness);


%%%output all data in form of tables and figures

%create a table of the essencial data
Parameter={'Major axis length'; 'Minor axis length';'Taylor parameter';'Eccentricity';'Flatness'};
Maximum=[maxMajorAxisLength;maxMinorAxisLength;maxTaylor;maxEccentricity;maxFlatness];
Minimum=[minMajorAxisLength;minMinorAxisLength;minTaylor;minEccentricity;minFlatness];
Max_frame=[I_maxMajorAxisLength;I_maxMinorAxisLength;I_maxTaylor;I_MaxEccentricity;I_maxFlatness];
Min_frame=[I_minMajorAxisLength;I_minMinorAxisLength;I_minTaylor;I_MinEccentricity;I_minFlatness];

T = table(Maximum,Minimum,Max_frame,Min_frame,...
    'RowNames',Parameter)

% plotting all info as a function of frame number. 
%% plotting major axis length

%------major axis---------------------------
FigMajorAxes= figure;
%create axis
axes1=axes('parent',FigMajorAxes,'FontSize',14);
box(axes1,'on');
hold(axes1,'on');

% Create xlabel
xlabel({'frame number'});
% Create ylabel
ylabel({'major axis length(pixel)'});
plot(LongAxis_value);
%----------------------------------------

%% plotting minor axis length
%------minor axis---------------------------
FigMinorAxes= figure;
%create axis
axes2=axes('parent',FigMinorAxes,'FontSize',14);
box(axes2,'on');
hold(axes2,'on');

% Create xlabel
xlabel({'frame number'});
% Create ylabel
ylabel({'minor axis length(pixel)'});
plot(shortAxis_value);
%----------------------------------------

%% plotting Orientation
%------Orientation---------------------------
FigOrientation= figure;
%create axis
axes3=axes('parent',FigOrientation,'FontSize',14);
box(axes3,'on');
hold(axes3,'on');

% Create xlabel
xlabel({'frame number'});
% Create ylabel
ylabel({'orientation'});
plot(orientation_value);
%----------------------------------------

%% plotting centroid location
%------centroid---------------------------
Figcentroid= figure;
%create axis
axes4=axes('parent',Figcentroid,'FontSize',14);
box(axes4,'on');
hold(axes4,'on');

% Create xlabel
xlabel({'frame number'});
% Create ylabel
ylabel({'centroid location (pixel)'});
plot1=plot(centroid_value_x);
set(plot1,'DisplayName','x');
hold on;
plot2=plot(centroid_value_y);
set(plot2,'DisplayName','y');
% Create legend
legend1 = legend(axes4,'show');
set(legend1,'FontSize',12.6);

%----------------------------------------

%% plotting Perimeter
%------Perimeter---------------------------
FigPerimeter= figure;
%create axis
axes5=axes('parent',FigPerimeter,'FontSize',14);
box(axes5,'on');
hold(axes5,'on');

% Create xlabel
xlabel({'frame number'});
% Create ylabel
ylabel({'Perimeter(pixel)'});
plot(Perimeter_value);
%----------------------------------------

%% plotting area
%------Area---------------------------
FigTaylorParam= figure;
%create axis
axes6=axes('parent',FigTaylorParam,'FontSize',14);
box(axes6,'on');
hold(axes6,'on');

% Create xlabel
xlabel({'frame number'});
% Create ylabel
ylabel({'Area(pixel)'});
plot(Area_value);

%% plotting Taylor Parameter

%------Taylor Parameter---------------------------
FigTaylorParam= figure;
%create axis
axes7=axes('parent',FigTaylorParam,'FontSize',14);
box(axes7,'on');
hold(axes7,'on');

% Create xlabel
xlabel({'frame number'});
% Create ylabel
ylabel({'Taylor Parameter'});
plot(TaylorParam);
