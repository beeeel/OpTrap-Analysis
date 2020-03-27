function [Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num)
%% [Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num)
% Load an Imstack, and according meta and info structs

validateattributes(Num,{'string','char'},{'nonempty','scalartext'},'Num')

[~, HName] = system('hostname');
HName = strsplit(HName);
if strcmp(HName{1}, 'will-linux')
    datadir = '/home/will/Documents/data/OpTrap/';
elseif strcmp(HName{1},'bowtie')
    datadir = '/home/ppxwh2/Documents/data/OpTrap';
end

infosdir = [datadir 'infos/'];
if isempty(whos('Imstack')); Imstack = {{0,''}}; end % Create an empty

if strcmp(CellType,'hela')
    S = load([infosdir 'info-hela_ctrl_s_020_tr_70_' Num '.mat']);
    imfile = [datadir '0610/Deformation/hela_ctrl_s_020_tr_70_' Num '_MMStack.ome.tif'];
    if ~strcmp(imfile, Imstack{1}{1,2}(1:length(imfile)))
        Imstack = load_imstack('0610/Deformation','ctrl','020','70',0,Num);
    end
elseif strcmp(Set,'Jenna')
    error('you need to update this branch of the load operation')
    S = load(strcat('/home/ppxwh2/Documents/data/OpTrap/infos/info_Jenna_test_no_erode',Num,'.mat'));
    imfile = ['/home/ppxwh2/Documents/data/OpTrap/1119_Jenna/191119_thp1_ctrl_s_010_tr_50_',Num,'_MMStack.ome.tif'];
    if ~strcmp(imfile, Imstack{1}{1,2}(1:length(imfile)))
        Imstack = bfopen(imfile);
    end
elseif strcmp(CellType,'HL60')
    if strcmp(Set,'with_drugs')
        matName = ['_190717_HL60_' Num '_0.020mm-1_1'];
    elseif strcmp(Set,'normoxia')
        matName = ['_100717_' Num '_' CellType '_1']; 
    end
    S = load(['/home/ppxwh2/Documents/data/OpTrap/infos/info_' CellType '_' Set matName '.mat']);
    imfile = ['/home/ppxwh2/Documents/data/OpTrap/2017_10_movies-from-aishah/'...
        CellType '/' CellType '_' Set '/' matName(2:end) '.avi'];
    if ~strcmp(imfile, Imstack{1}{1,2}(1:length(imfile)))
        Imstack = avi_to_imstack(imfile);
    end
elseif strcmp(CellType,'LS174T')
    imfile = ['200717_' Num '_LS174T_' Set '_1.avi'];
    if strcmp(Set,'hypoxia'); imfile(2) = '1'; end
    S = load([infosdir 'info_seg_LS174T_' Set '_' imfile(1:end-4) '.mat']);
    Imstack = avi_to_imstack([datadir '2017_10_movies-from-aishah/LS174T/' imfile]);
end
info = S.info;
meta = S.meta;
end