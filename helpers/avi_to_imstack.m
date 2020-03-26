function Imstack = avi_to_imstack(filename)
% Load an avi file, return a cell array formatted for use with
% PostProcessCellDeform_v2 etc
% Imstack = avi_to_imstack(filename)
% filename must be absolute path, or in a directory on MATLAB's path

v = VideoReader(filename);
vid = read(v, [1 inf]);
N_frames = size(vid, 4);
file_path = [v.Path '/' v.Name];

Imstack = cell(1);
Imstack{1} = cell(N_frames,2);

for frame = 1:N_frames
    if strcmp(v.VideoFormat, 'RGB24')
        Imstack{1}{frame,1} = rgb2gray(squeeze(vid(:,:,:,frame)));
    elseif strcmp(v.VideoFormat,'Grayscale')
        Imstack{1}{frame,1} = squeeze(vid(:,:,:,frame));
    end
    Imstack{1}{frame,2} = [file_path ' ' num2str((frame - 1) / v.FrameRate) 's'];
    ProgressBar(frame/N_frames)
end

end