function compare_info_meta_imstack(Info, Meta, Imstack)
SplitInfo = strsplit(Info,'/');
SplitMeta = strsplit(Meta, '/');
SplitImstack = strsplit(Imstack,{'/',' '});

if isempty(whos('Imstack'))
    error('You need an Imstack')
elseif isempty(whos('info')) || isempty(whos('meta'))
    error('You need info and meta')
elseif ~strcmp(SplitInfo(1).filepath,SplitMeta.filepath) || ~strncmp(SplitMeta.filepath, SplitImstack{1}{1,2}, length(SplitMeta.filepath))
    fprintf('Info for dataset %s\n Meta for dataset %s\n Imstack for dataset %s\n',...
        SplitInfo(1).filepath, SplitMeta.filepath, SplitImstack{1}{1,2}(1:length(SplitMeta.filepath)))
    error('info, meta, and Imstack need to be for the same dataset')
end
end
