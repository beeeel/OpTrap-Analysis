function compare_info_meta_imstack(Info, Meta, Imstack)
%% Check an info and meta are from the same Imstack as yours

% The filename is in the final element of info and meta, but Imstack can
% has a time for the frame after 
SplitInfo = strsplit(Info(1).filepath,'/');
SplitMeta = strsplit(Meta.filepath, '/');

SplitImstack = strsplit(Imstack{1}{1,2},{'/',' '});

if length(SplitImstack{end}) == 2
    SplitImstack = SplitImstack(1:end-1);
end

if isempty(whos('Imstack'))
    error('You need an Imstack')
elseif isempty(whos('Info')) || isempty(whos('Meta'))
    error('You need info and meta')
elseif ~strcmp(SplitInfo{end},SplitMeta{end}) || ~strcmp(SplitMeta{end}, SplitImstack{end})
    fprintf('Info for dataset %s\n Meta for dataset %s\n Imstack for dataset %s\n',...
        SplitInfo(1).filepath, SplitMeta.filepath, SplitImstack{1}{1,2}(1:length(SplitMeta.filepath)))
    error('info, meta, and Imstack need to be for the same dataset')
end
end
