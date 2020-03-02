function PostProcessCellDeformBatch_v2(dirIn, saveall)
%       PostProcessCellDeformBatch_v2(dirIn, saveall)
%% Batch load, process, and save of cell data
% Sequentially load, process, and save the output of cell imaging data. If
% save flag is 1, save the whole info structure, otherwise just save a
% summary figure (major + minor axis lengths, area and find_cell radius),
%
%
% This will overwrite any files with the same save name 
% Directory is chosen with dirIn, and all files in that directory will be
% iterated over. Only files ending with '.tif' will be loaded
% Output directory is OpTrap/infos, plots are in OpTrap/processing_plots.

% Parameters for PPCD
find_cell = {'Lsigma', 500, 'Lalpha', 0.95 };
segment_cell = {'ellipseFitVal',0,'iterations',250};
seg_cell_v = 5;
find_cell_v = 2;

% Common directories for all treatments:
outpath = '/home/ppxwh2/Documents/data/OpTrap/infos/';
rootdir = '/home/ppxwh2/Documents/data/OpTrap/';
% Define the directories we're looking at, and get the list of files in our
% raw data directory
dir = strcat(rootdir, dirIn);

% Later bits need the directory to end with a /
if ~strcmp(dir(end), '/')
    dir = strcat(dir, '/');
end

filelist = strsplit(ls(dir));
% Iterate over the list of files
for file = filelist(1:end-1)
    SplitFile = strsplit(file{:}, '.');
    if sum(strcmp(SplitFile{end}, {'tif','avi'})) == 1
        % Load in the raw data file
        disp(strcat('Loading file: ', dir, file))
        switch SplitFile{end}
            case 'tif'
                Imstack = bfopen(strcat(dir, file{:}));
            case 'avi'
                Imstack = avi_to_imstack([dir file{:}]);
        end
        % Process it with PPCD
        info = PostProcessCellDeform_v2(Imstack, 'seg_cell_v',seg_cell_v,...
            'find_cell_v',find_cell_v,'segment_cell',segment_cell)
        
        if saveall == 1
            % Create a filename to save - filenames end in MMStack.ome.tif (16char)
            % which is stripped by the indexing
            fileOut = strcat(outpath,'info-otsu_', file{:}(1:end-16))
            % If the struct is larger than 2GB, it will need to be saved with v7.3,
            % but I've had some bugs loading v7.3 mat files, so don't save it like
            % that if I don't have to.
            % Get properties for the info variable
            vars = whos('info');
            % If the size of the info struct is bigger than 1.9GB, use v7.3
            if vars.bytes > 1.9e9
                save_ver = '-v7.3';
            else
                save_ver = '-v7';
            end
            saveall(fileOut, 'info', save_ver)
        end
        
        % Make a basic set of 4 plots and save the file
        figure(1)
        subplot(221)
        plot([info.MajorAxisLength])
        title('Major Axis')
        subplot(222)
        plot([info.MinorAxisLength])
        title('Minor Axis')
        subplot(223)
        plot([info.Area])
        title('Area')
        subplot(224)
        plot([info.radius])
        title('Radius (find_cell)')
        filename = strsplit(info(1).filepath,{'/','.'});
        hgsave(strcat('~/Documents/data/OpTrap/processing_plots/',...
            filename{end-2},'_summaryplot'))
    else
        disp(strcat('Skipping file: ', dir, file))
    end
end

end
