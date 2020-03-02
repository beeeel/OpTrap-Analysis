function [mask] = ellipse_mask(semi_x,semi_y,centre_x,centre_y,mask)
mask_size = size(mask);

%(x./xa)^2 + (y/ya)^2 = 1
mask(:) = 0;

for i = 1:mask_size(1)
    for j = 1:mask_size(2)
        if ((i - centre_y)./semi_y).^2 + ((j - centre_x)./semi_x).^2 <= 1
            mask(i,j) = 1;
        end
    end
end
end