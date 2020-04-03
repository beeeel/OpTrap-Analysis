function [Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num, varargin)
%% [Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num)
% Load an Imstack, and according meta and info structs

validateattributes(Num,{'string','char'},{'nonempty','scalartext'},'Num')

if nargin == 3
    
    
end

[~, HName] = system('hostname');
HName = strsplit(HName);
if strcmp(HName{1}, 'will-linux')
    DataDir = '/home/will/Documents/data/OpTrap/';
elseif strcmp(HName{1},'bowtie')
    DataDir = '/home/ppxwh2/Documents/data/OpTrap/';
end

InfosDir = [DataDir 'infos/'];
if isempty(whos('Imstack')); Imstack = {{0,''}}; end % Create an empty

if strcmp(CellType,'hela')
    S = load([InfosDir 'info-hela_ctrl_s_020_tr_70_' Num '.mat']);
    Imstack = load_imstack('0610/Deformation','ctrl','020','70',0,Num);
elseif strcmp(Set,'Jenna')
    S = load(strcat(InfosDir, 'info_Jenna_test_no_erode',Num,'.mat'));
    imfile = [DataDir '1119_Jenna/191119_thp1_ctrl_s_010_tr_50_',Num,'_MMStack.ome.tif'];
    Imstack = bfopen(imfile);
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
    S = load([InfosDir 'info_reduced_seg_' CellType '_' Set matName '.mat']);
    imfile = [DataDir '2017_10_movies-from-aishah/'...
        CellType '/' CellType '_' Set '/' matName(2:end) '.avi'];
    Imstack = avi_to_imstack(imfile);
elseif strcmp(CellType,'LS174T')
    imfile = ['200717_' Num '_LS174T_' Set '_1.avi'];
    if strcmp(Set,'hypoxia'); imfile(2) = '1'; end
    S = load([InfosDir 'info_seg_LS174T_' Set '_' imfile(1:end-4) '.mat']);
    Imstack = avi_to_imstack([DataDir '2017_10_movies-from-aishah/LS174T/' imfile]);
end
info = S.info;
meta = S.meta;
end