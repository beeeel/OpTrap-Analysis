function view_stack( Imstack, varargin)
%view_stack( Imstack, [crop], [pause], [repeats])
%   Run through a stack of images, showing a box with edges defined by
%   crop. Time between frames defined in sec by pause, number of cycles
%   through whole stack defined by repeats. 
%   Default parameters are no crop box, 0.1s pause, 1 repeat.

n_frames = size(Imstack{1},1);

% defaults
p_t = 0.5;
repeats = 1;
crop = [];

if nargin > 1
        % Crop will be part of the cell_data struct outputted from
        % find_cell
        if isstruct(varargin{1})
            crop = zeros(4, n_frames);
            for frame = 1:n_frames
                crop(:, frame) = varargin{1}(frame).crop;
            end
        elseif sum(size(varargin{1}) == [4, 1]) == 2
            crop = zeros(4, n_frames);
            for frame = 1:n_frames
                crop(:,frame) = varargin{1};
            end
        elseif sum(size(varargin{1}) == [4, n_frames]) == 2
            crop = varargin{1};
        end
        % Pause time will be a number which is not an integer
        if sum(size(varargin{2}) == [1 1]) == 2
            if isnumeric(varargin{2}) == 1 && ...
                    round(varargin{2}) ~= varargin{2};
                p_t = varargin{2};
            end
        end
        % Repeats will be an integer
        if isinteger(varargin{3}) == 1;
            repeats = varargin{3};
        end
end
% Show values for checking
disp(['Size of crop array = ', size(crop)])
disp(['Number of repeats = ', repeats])
disp(['Pause time = ', p_t, 's'])


figure
for rep = 1:repeats
    for frame = 1:n_frames
        % For each frame, display it, draw the crop box, write the title
        imagesc(Imstack{1}{frame,1})
        axis image off, colormap gray
        if size(crop,1) > 1
            hold on
            plot([crop(1, frame), crop(2, frame), crop(2, frame), crop(1, frame), crop(1, frame)], ...
                [crop(3, frame), crop(3, frame), crop(4, frame), crop(4, frame), crop(3, frame)], ...
                'r--')
        end
        title(['Frame ', num2str(frame)])
        pause(p_t);
    end
end
end

