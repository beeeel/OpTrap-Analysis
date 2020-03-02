outpath = '/home/ppxwh2/Documents/data/OpTrap/infos/';

for rep = ['3','5','6']
    file = strcat('hela_ctr_s_020_tr_70_',rep);
    if ~strcmp(rep,'3')
        Imstack = load_imstack('0610/Deformation','ctrl','020','70',0,'6');
    
    info = PostProcessCellDeform_v2(Imstack,'seg_cell_v',3,'find_cell_v',2);
    % Create a filename to save - filenames end in MMStack.ome.tif (16char)
    % which is stripped by the indexing
    fileOut = strcat(outpath,'info_', file)
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
    save(fileOut, 'info', save_ver)
    end
end