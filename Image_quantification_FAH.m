%{
Project name: Image Processing for Tracking an Optically Trapped Cell
Programmer  : Felix A.Halim

The following code calculates surface drift, cell deformation and motion
Binary converstion is utilised for tracking surface drift
Edge detection is utilised to indicate deformation
The pattern matching is utilised to quantify motion
%}

close all
clear variables
clc

%% variables input 

frames=130;         % number of frames
px_to_mc=10;        % resoliton (pixel per micron)
filename = '~/Documents/data/OpTrap/0610/Deformation/hela_ctrl_s_020_tr_70_1_MMStack.ome';     % video names
offset=0;           % surface offset (boost) in pixel
Type=1;             % surface elimination type (1 or 2)

%% advance input

%smoothing setting
sigma=0.2;          % 0 to 1 amplitude of edges
alpha=5;            % higher more smooth

%pattern mathing setting
[optimizer,metric]= imregconfig('monomodal'); %setting for patern matching

%video output setting
%v=VideoWriter(strcat(filename,' animation.avi'),'Motion JPEG AVI');

%graph position
graph_pos=[680 0 560 1000]; %this setting is based on my screen (1920x1080)

%% ROI selection

%open (v); %open the address
a=1;
I_O=imread(strcat(filename,'.tif'),a);  %read the first image
I_O=imreducehaze(I_O);                  %image haze removal
figure,imshow(I_O), title('Original Image')
mask=drawrectangle;                     %interactive action of cropping
crop=createMask(mask);                  %create mask
[B,L,N]=bwboundaries(crop,'noholes');   %identify the mask
property=regionprops(L,'BoundingBox');  %extract outline/position of the mask
ROI=property.BoundingBox;               
Size=[ROI(3) ROI(4)];                   %save the size of cropped image
parameter=max(Size);                    %save as threshold for later

%% Preallocate memory
BBoxL = zeros(frame,1);
BBoxR = zeros(frame,1);
BBoxT = zeros(frame,1);
BBoxB = zeros(frame,1);
Surf_H = zeros(frame,1);
Area = zeros(frame,1);
simx = zeros(frames - 1,1);
simy = zeros(frames - 1,1);
simb = zeros(frames - 1,1);
simc = zeros(frames - 1,1);

%% data collection
for a=1:frames                          %loop for all frames
%% image preprocess
    
I_O=imread(strcat(filename,'.tif'),a);   %read file
I_O=imcrop(I_O,ROI);                     %take the crop ROI
I_O=imreducehaze(I_O);                   %image haze removal
I_Ori=I_O;                               %back up image

I_S=locallapfilt(I_Ori,sigma, alpha);    %smoothed iamge 
                                         %I_S is used on Edge Detection
%% filtering and morphological process section

% surface binary
I_O=imbinarize(I_O);                %binary conversion
I_O=bwareaopen(I_O,parameter);      %remove small noises (i.e grain)
se=strel('disk',1);                 %patch creation               
I_O=imclose(I_O,se);                %fill the holes
I_O=imfill(I_O,'holes');            %morpological process
                                    %I_O is clean binary image
%% Identification/ Acquisation

xx=1;                                     %tracking dummy
[B1,L1,N1]=bwboundaries(1-I_O,'noholes'); %identify each detected objects
property1=regionprops(L1,'all');          %extract all available features from the objects

% check the biggest Area (surface object)
Area(a)=property1(1).Area;              %put the area starter (first) 
if N1>1                                 %if there is more than one object 
    for na=2:N1
        if property1(na).Area>Area(a)   %search for the larger area
        Area(a)=property1(na).Area;     %update the new bigges Area
        xx=na;                          %update the tracker
        end              
    end   
end

% vertical or horizontal (orientation and location)
BB1=property1(xx).BoundingBox;          %record the bounding box
I_FO1=1-I_O;                            %I_FO1 is the first method of surface mask
I_FO2=1-I_O;                            %I_FO2 is the second method of surface mask
if property1(xx).BoundingBox(3)==Size(1)+1  %if its horizontal
   hv=4;        %indicator of vertical
        % Up or Down
        if property1(xx).Centroid(2)>Size(2)/2            
            sprintf('bottom');surf_loc=1;
            I_FO1(ceil(BB1(2)-offset):ceil((BB1(2)+BB1(4)-offset-1)),ceil(BB1(1)):ceil(BB1(1)+BB1(3)-1))=property1(xx).FilledImage;
            I_FO1(ceil((BB1(2)+BB1(4)-offset-1)):Size(2)+1,:)=ones;
            I_FO2(ceil(BB1(2)-offset):ceil((BB1(2)+BB1(4)-offset-1)),ceil(BB1(1)):ceil(BB1(1)+BB1(3)-1))=ones;
        else
            sprintf('top');surf_loc=2;
            I_FO1(ceil(BB1(2)+offset):ceil((BB1(2)+BB1(4)+offset-1)),ceil(BB1(1)):ceil(BB1(1)+BB1(3)-1))=property1(xx).FilledImage;
            I_FO1(1:ceil(BB1(2)+offset),:)=ones;
            I_FO2(ceil(BB1(2)+offset):ceil((BB1(2)+BB1(4)+offset-1)),ceil(BB1(1)):ceil(BB1(1)+BB1(3)-1))=ones;
        end
else                                          %if vertical
   hv=3;       %indicator of horizontal
        % Left or Right
        if property1(xx).Centroid(1)<Size(1)/2
            sprintf('left');surf_loc=3;
            I_FO1(ceil(BB1(2)):ceil(BB1(2)+BB1(4)-1),ceil(BB1(1)+offset):ceil(BB1(1)+BB1(3)+offset-1))=property1(xx).FilledImage;
            I_FO1(:,1:ceil(BB1(1)+offset))=ones;
            I_FO2(ceil(BB1(2)):ceil(BB1(2)+BB1(4)-1),ceil(BB1(1)+offset):ceil(BB1(1)+BB1(3)+offset-1))=ones;
        else
            sprintf('right');surf_loc=4;
            I_FO1(ceil(BB1(2)):ceil(BB1(2)+BB1(4)-1),ceil(BB1(1)-offset):ceil(BB1(1)+BB1(3)-offset-1))=property1(xx).FilledImage;
            I_FO1(:,ceil(BB1(1)+BB1(3)-offset-1):Size(1)+1)=ones;
            I_FO2(ceil(BB1(2)):ceil(BB1(2)+BB1(4)-1),ceil(BB1(1)-offset):ceil(BB1(1)+BB1(3)-offset-1))=ones;
        end
end

%surface drift data saving
switch surf_loc
    case 1    %case the surf is at bottom
        Surf_H(a)=BB1(2);
    case 2    %case the surf is at top
        Surf_H(a)=BB1(2)+BB1(4);
    case 3    %case the surf is at L
        Surf_H(a)=BB1(1)+BB1(3);
    case 4    %case the surf is at R
        Surf_H(a)=BB1(1);
    otherwise
        sprintf ('unk surf_lock')
end

%surface mask selection
if Type==2      
    I_F=I_FO2;  %second method surface mask 
else
    I_F=I_FO1;  %first method surface mask
end

%end of surface data procesisng and sacing
%% Cell Acquisation

% smoothed image
I_S=edge(I_S,'canny');                %apply canny edge detection
I_S=immultiply(1-I_F,I_S);            %exclude the surface

% Boundary from smoothed Image
[SB,SL,SN]=bwboundaries(I_S,'noholes'); %identify each detected objects/lines
propertyS=regionprops(SL,'all');        %extract all available features from the objects/lines

%put the lines on the edges of cropped images
BBoxLT=Size(1);
BBoxRT=0;
BBoxTT=Size(2);
BBoxBT=0;

% obtaining the boundary
for kS=1:length(SB)
    if propertyS(kS).Area > 10         %if the lenght is long enought to be processed
        BBox=propertyS(kS).BoundingBox;     %record the BOunding box    
        
        %get the left boundary
        if BBox(1)<BBoxLT && BBox(1)>ROI(4)/10
        BBoxLT=BBox(1);
        end
        
        %get the right boundary
        BBoxTemp=BBox(1)+BBox(3);
        if BBoxTemp>BBoxRT && (BBoxTemp)<(Size(1)-ROI(4)/10)
        BBoxRT=BBoxTemp;
        end
        
        % get the top boundary
        if BBoxTT>BBox(2) && BBox(2)>ROI(3)/10
        BBoxTT=BBox(2);
        end
        
        %get the bottom boundary
        BBoxTemp2=BBox(2)+BBox(4);
        if BBoxBT<BBoxTemp2 && BBoxTemp2<(Size(2)-ROI(3)/10)
        BBoxBT=BBox(2)+BBox(4);
        end
    end
end

%% save values

%assign the borders to another names (global names)
BBoxL(a)=BBoxLT;
BBoxR(a)=BBoxRT;
BBoxT(a)=BBoxTT;
BBoxB(a)=BBoxBT;

% plotting images
annotation = sprintf('time %i s',a*6);      %time indicator (top left of vdieo)
I_Ori=insertText(I_Ori,[1 1],annotation);
I_Ori=insertShape(I_Ori,'Rectangle',[BBoxLT BBoxTT BBoxRT-BBoxLT BBoxBT-BBoxTT]); %cell box drawing
imshow(I_Ori)

%M=getframe(gcf);    %put tup the frame to be recorded
%writeVideo(v,M)     %video recording each M frame
 
end

%close(v)    %end of recording
%% global data saving 

x=1:frames;
yBBW=(BBoxR(x)-BBoxL(x))/px_to_mc;
yBBT=(BBoxB(x)-BBoxT(x))/px_to_mc;
ySH=Surf_H/px_to_mc;

%ratio calculation
if(yBBW>=0)
yR(x)=yBBT(x)./yBBW(x);
else
yR(x)=0;    
end

%% pattern changes

%get the dimension for the cell region
W=max(BBoxR)-min(BBoxL);
H=max(BBoxB)-min(BBoxT);
ROI_plus=[min(BBoxL) min(BBoxT) W H];
ROI_T=[ROI(1) ROI(2) 0 0]+ROI_plus;

for a=1:frames-1
   
I_O1=imread(strcat(filename,'.tif'),a);     %read file
I_O1=imcrop(I_O1,ROI_T);                    %take the ROI
I_O1=imreducehaze(I_O1);
I_Ofix=I_O1;                                %back up image

I_O2=imread(strcat(filename,'.tif'),a+1);   %read file
I_O2=imcrop(I_O2,ROI_T);                    %take the ROI
I_O2=imreducehaze(I_O2);
I_Omov=I_O2;                                %back up image

%% pattern mathinh 
tform2 = imregtform(I_Omov, I_Ofix, 'rigid', optimizer, metric);
simx(a)=tform2.T(3,1);  %x translation
simy(a)=tform2.T(3,2);  %y translation
simb(a)=tform2.T(1,2);  %'b' value: rotation angle 1
simc(a)=tform2.T(2,1);  %'c' value: rotation angle 2

end

%% global data saving

x2=1:frames-1;

yxs=simx(x2);
yys=simy(x2);
yrb=asind(simb(x2));

%% data plotting

f=figure;
ax7=subplot(7,1,1);
plot(x*6,yBBW,'b')
title(ax7,'Cell Width')
xlabel(ax7,'time(s)') 
ylabel (ax7,'widht(um)')

ax8=subplot(7,1,2);
plot(x*6,yBBT,'r')
title(ax8,'Cell Height')
xlabel(ax8,'time(s)') 
ylabel (ax8,'height(um)')

ax5=subplot(7,1,3);
plot(x*6,yR,'k')
title(ax5,'H/W Ratio')
xlabel(ax5,'time(s)') 
ylabel (ax5,'ratio')

ax6=subplot(7,1,4);
plot(x*6,ySH,'g')
title(ax6,'Surface Drift')
xlabel(ax6,'time(s)') 
ylabel (ax6,'position(um)')

ax2=subplot(7,1,5);
plot(x2*6,yxs,'b')
title(ax2,'x translation')
xlabel(ax2,'time(s)') 
axis([1 (frames-1)*6 -4 4])
ylabel (ax2,'changes(um)')

ax4=subplot(7,1,6);
plot(x2*6,yys,'r')
title(ax4,'y translation')
xlabel(ax4,'time(s)') 
axis([1 (frames-1)*6 -4 4])
ylabel (ax4,'changes(um)')

ax6=subplot(7,1,7);
plot(x2*6,yrb,'g')
title(ax6,'rotation')
xlabel(ax6,'time(s)') 
axis([1 (frames-1)*6 -2 2])
ylabel (ax6,'changes (degree)')

f.Position=graph_pos;
%f.WindowState='maximize'; %enable if the graph wanted to be in fullscreen
%saveas(gcf,strcat(filename,'_result','.jpg'))