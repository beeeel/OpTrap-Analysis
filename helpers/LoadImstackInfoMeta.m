function LoadImstackInfoMeta(CellType, Set, Num, varargin)
%% [Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num)
% Load an Imstack, and according meta and info structs

validateattributes(Num,{'string','char'},{'nonempty','scalartext'},'Num')

global Imstack info meta
persistent InfoFileOld DataFileOld

[~, HName] = system('hostname');
HName = strsplit(HName);
if strcmp(HName{1}, 'will-linux')
    DataDir = '/home/will/Documents/data/OpTrap/';
elseif strcmp(HName{1},'bowtie')
    DataDir = '/home/ppxwh2/Documents/data/OpTrap/';
end

InfosDir = [DataDir 'infos/'];

% 1) Get filenames for data and info
if strcmp(CellType,'hela')
    InfoFile = {[InfosDir 'info-hela_ctrl_s_020_tr_70_' Num '.mat']};
    DataFile = {'0610/Deformation','ctrl','020','70',0,Num};
    Loader = 'load_imstack';
elseif strcmp(Set,'Jenna')
    InfoFile = {strcat(InfosDir,'info_Jenna_test_no_erode',Num,'.mat')};
    DataFile = {[DataDir '1119_Jenna/191119_thp1_ctrl_s_010_tr_50_',Num,'_MMStack.ome.tif']};
    Loader = 'bfopen';
elseif strcmp(CellType,'HL60')
    if strcmp(Set,'with_drugs')
        matName = ['_190717_HL60_' Num '_0.020mm-1_1'];
    elseif strcmp(Set,'normoxia')
        if str2double(Num) < 11
            matName = ['_100717_' Num '_' CellType '_1'];
        else
            matName = ['_' CellType '_' Num '_0.020mms-1_1'];
        end 
    end
    InfoFile = {[InfosDir 'info_reduced_seg_' CellType '_' Set matName '.mat']};
    DataFile = {[DataDir '2017_10_movies-from-aishah/'...
        CellType '/' CellType '_' Set '/' matName(2:end) '.avi']};
elseif strcmp(CellType,'LS174T')
    switch Set
        case 'hypoxia'
            Date = '210717';
        case 'normoxia'
            Date = '200717';
    end
    DataFile = {[DataDir '2017_10_movies-from-aishah/LS174T/' Date '_' Num '_LS174T_' Set '_1.avi']};
    InfoFile = {[InfosDir 'info_seg_LS174T_' Set '_' Date '_' Num '_LS174T_' Set '_1.mat']};
end

% 2) Load it if the last file was different to this
if ~strcmp(InfoFileOld, InfoFile) || ~strcmp(DataFileOld, DataFile)
    Imstack = eval([Loader '(DataFile{:});']);
    S = load(InfoFile{:});
    % 3) Extract info and meta for output
    info = S.info;
    meta = S.meta;
    % 4) Update trackers
    InfoFileOld = InfoFile;
    DataFileOld = DataFile;
end
end