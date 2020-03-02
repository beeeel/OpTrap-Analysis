function [ r_out] = cell_radii( Imstack, crop)
%cell_radii measures cell radius in cropped image using bright halo
%   Measure radii from the bright halo
r_out = zeros(2, size(Imstack{1}, 1));
%warndlg('work in progress!')
n_fr = size(Imstack,1);
for frame = 1:n_fr
    if isnan(Imstack{frame})
        r_out(frame,:) = [nan, nan];
    else
        w = crop() - crop();
        h = crop() - crop();
        
        r_x = 0;
        r_y = 0;
        
        for Hslice = 1:w
            % Get the indices where the value is the max for the first and second
            % halves of the image
            
            % Currently indexes half the width of the image - need to index for an
            % area taken using the radius and centre of the cirle
            first_peak = find(Imstack{frame}(Hslice,:) == ...
                max(Imstack{frame}(Hslice,ceil(1:w/2))));
            second_peak = find(Imstack{frame}(Hslice,:) == ...
                max(Imstack{frame}(Hslice,ceil(w/2:end))));
            if second_peak(1) - first_peak(1) > r_x
                r_x = second_peak(1) - first_peak(1);
            end
        end
        
        for Vslice = 1:h
            % Get the indices where the value is the max for the first and second
            % halves of the image
            
            % Currently indexes half the height of the image - need to index for an
            % area taken using the radius and centre of the cirle
            first_peak = find(Imstack{frame}(:, Vslice) == ...
                max(Imstack{frame}(ceil(1 :w/2),Vslice)));
            second_peak = find(Imstack{frame}(:,Vslice) == ...
                max(Imstack{frame}(ceil(w/2:end),Vslice)));
            if second_peak(1) - first_peak(1) > r_y
                r_y = second_peak(1) - first_peak(1);
            end
        end
        r_out(frame, :) = [r_x, r_y];
    end
end

end

