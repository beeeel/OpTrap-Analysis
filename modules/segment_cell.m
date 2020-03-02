function [mask] = segment_cell(im,varargin)
% SEGMENT_CELL  Takes an input grayscale image and generates a binary image mask.
%   Mask = SEGMENT_CELL(IM, OPTIONAL:MinSize, OPTIONAL:EllipseFitVal).
%
%   Segments an input grayscale image given the image and minimum size
%   of cell, returns a Binary image mask. If Ellipse Fit Val is set equal
%   to 1 then Josh's Ellipse fitting function will be used to fit an
%   ellipse to the mask, otherwise just the mask will be returned. 
%   
%
%   By Default, EllipseFitVal = 0
%

% default values
if nargin == 1
    minSize = 12000; % Number of pixels
    ellipseFitVal = 0;
elseif nargin == 2
    minSize = varargin{1};
    ellipseFitVal = 0;
else
    minSize = varargin{1};
    ellipseFitVal = varargin{2};
    if minSize < ellipseFitVal
        error('Input arguments incorrect order');
    end
end



%first check that Image is grayscale, commented out for now
%only works for grayscal image
%if 
    %Finding threshold level, fudging this and thresholding with fudge
    [~, threshold] = edge(im, 'canny');
    fudgeFactor = .5;
    bwIm = edge(im,'canny', threshold * fudgeFactor);
    
    %Structuring elements for dilation
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    
    %Dilate image with line structuring elements
    bwImDilate = imdilate(bwIm, [se90 se0]);
    
    %fill holes in objects
    bwImFilled = imfill(bwImDilate, 'holes');
    
    %get rid of objects touching the border as the cell will be in center of
    %image
    bwImNoBorders = imclearborder(bwImFilled, 4);
    
    %Erode edges with diamond structuring element
    seD = strel('diamond',1);
    bwImEroded = imerode(bwImNoBorders,seD);
    bwImEroded = imerode(bwImEroded,seD);
    
    %Remove objects smaller than certain size pixels
    bwImNoSmallObj = bwareaopen(bwImEroded, minSize);
    
    
    %Now apply a kernal with window sized 51
    windowSize = 51;
    kernel = ones(windowSize) / windowSize ^ 2;
    bwImNoSmallObj = conv2(single(bwImNoSmallObj), kernel, 'same');
    bwSmoothed = bwImNoSmallObj > 0.5; % Rethreshold
    
    if ellipseFitVal == 1
        %Now generate Elipse paramaters and fit using Elipse fitting function
        [semiY,semiX,centreX,centreY] = ellipse_fit_fun_final(bwSmoothed);
        
        mask = ellipse_mask(semiX,semiY,centreX,centreY,bwSmoothed);
        
    elseif ellipseFitVal ==2
        %Else if ellipseFitVal ==2, use ellipseDetection method
        [bestFits] = ellipseDetection(bwIm);
        
        mask = ellipse_mask(majorAxis,minorAxis,centreX,centreY,bwSmoothed);
    else
        %If not fitting Ellipse, just set to previous Mask (BWSmoothed)
        mask = bwSmoothed;
    end
    mask = uint8(mask);
end
