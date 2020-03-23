function compare_info_meta_imstack(info, meta, Imstack)
if isempty(whos('Imstack'))
    error('You need an Imstack')
elseif isempty(whos('info')) || isempty(whos('meta'))
    error('You need info and meta')
elseif ~strcmp(info(1).filepath,meta.filepath) || ~strncmp(meta.filepath, Imstack{1}{1,2}, length(meta.filepath))
    fprintf('Info for dataset %s\n Meta for dataset %s\n Imstack for dataset %s\n',...
        info(1).filepath, meta.filepath, Imstack{1}{1,2}(1:length(meta.filepath)))
    error('info, meta, and Imstack need to be for the same dataset')
end
end
