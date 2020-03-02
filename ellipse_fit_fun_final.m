% Ellipse fitting function for fitting ellipses to binary segmentation
% masks.
% Author: Joshua McAteer | Date: 15/05/2018
% University of Nottingham,Faculty of Engineering,Optics and Photonics Group

function [semi_y,semi_x,centre_x_best,centre_y_best] = ellipse_fit_fun_final(mask)

%variables used for plotting and saving figures, set =1 to save/plot
plot_var = 0;
save_var = 0;

%Normalise input mask
mask = mask - min(mask(:));
mask = mask./max(mask(:));

mask_size = size(mask);

%project mask onto x and y axis by summing all pixels along axis
yproj = sum(mask,2);
xproj = sum(mask,1);

max_yproj = max(yproj);
min_yproj = min(yproj);

max_xproj = max(xproj);
min_xproj = min(xproj);

%threshold to cut out outer pixels
threshold = 0.3;

%set seperate thresholds for x and y axis, [max.*(1-thresh) + min]
ythreshold = max_yproj.*(1-threshold) + min_yproj;
xthreshold = max_xproj.*(1-threshold) + min_xproj;

%initilize range variables
yrange = [0,0];
xrange = [0,0];

%fill these variables with the max and min values to use after thresholding

%First find min by moving along y axis until pixel is greater than thresh
for i = 1:mask_size(1)
    if yproj(i) >=ythreshold
        yrange(1) = i;
        break
    end
end

%Flip the axis and move along reverse until pixel is greater than thresh
for i = fliplr(1:mask_size(1))
    if yproj(i) >=ythreshold
        yrange(2) = i;
        break
    end
end

for i = 1:mask_size(2)
    if xproj(i) >=xthreshold
        xrange(1) = i;
        break
    end
end
for i = fliplr(1:mask_size(2))
    if xproj(i) >=xthreshold
        xrange(2) = i;
        break
    end
end

%initial semi_major axis vals = range/2
semi_y = (yrange(2) - yrange(1))./2;
semi_x = (xrange(2) - xrange(1))./2;

%initial center vals, (max+min)/2
centre_y = (yrange(2) + yrange(1))./2;
centre_x = (xrange(2) + xrange(1))./2;

%set guess variables equal to the initial guess
semi_y_guess = semi_y;
semi_x_guess = semi_x;

centre_y_guess = centre_y;
centre_x_guess = centre_x;

%
change = 10;
a = mean(mask_size)./10;
c = mean(mask_size)./10;

a_decay = 0.99; % can measure precision by a*a_decay^1000
c_decay = 0.99;

%make initial guess and set initial goodness value
mask_guess = ellipse_mask(semi_x,semi_y,centre_x,centre_y,mask);

goodness = sum(abs(mask(:) - mask_guess(:)));

%make continous guesses for i =1:1000
%NOTE: could edit this to make more efficient, have some get out statement
%that runs if change is minimal etc.

for i = 1:1000
    mask_guess = ellipse_mask(semi_x_guess,semi_y_guess,centre_x_guess,centre_y_guess,mask);
    
    goodness_guess = sum(abs(mask(:) - mask_guess(:)));
    
    if goodness_guess < goodness
        
        goodness = goodness_guess;
        semi_x = abs(semi_x_guess);
        semi_y = abs(semi_y_guess);
        centre_y = abs(centre_y_guess);
        centre_x = abs(centre_x_guess); 
    end
    
    semi_x_guess = semi_x + a.*randn(1);
    semi_y_guess = semi_y + a.*randn(1);
    centre_x_guess = centre_x + c.*randn(1);
    centre_y_guess = centre_y + c.*randn(1);
    
    a = a.*a_decay;
    c = c.*c_decay;
    
end

%set final semi_major and semi minor
semi_major = semi_y;
semi_minor = semi_x;

centre_x_best = centre_x;
centre_y_best = centre_y;

if semi_x >= semi_y
    semi_major = semi_x;
    semi_minor = semi_y;
end

if save_var == 1
    saveas(gcf,['plot_',num2str(round((semi_major./semi_minor).*100)),'.png']);
end
end

