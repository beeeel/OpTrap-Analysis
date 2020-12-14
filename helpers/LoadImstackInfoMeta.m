function varargout = LoadImstackInfoMeta(CellType, Set, Num, varargin)
%% LoadImstackInfoMeta(CellType, Set, Num, [ImstackOnly, InfoPrefix, InfoSuffix, SimpleMatName])
% Load an Imstack, and (maybe) according meta and info structs into global
% variables

% -1) Check inputs are correct
validateattributes(Num,{'string','char'},{'nonempty','scalartext'},'Num')
if nargin > 3
    validateattributes(varargin{1},{'logical'},{'nonempty','scalar'},'ImstackOnly')
end
if nargin > 4
    validateattributes(varargin{2},{'string','char'},{'scalartext'},'InfoPrefix')
end
if nargin > 5
    validateattributes(varargin{3},{'string','char'},{'nonempty','scalartext'},'InfoSuffix')
end
if nargin > 6
    validateattributes(varargin{4},{'logical'},{'nonempty','scalar'},'SimpleMatName')
end

% 0) Prepare to get filenames
if nargout == 0
    global Imstack info meta
end
persistent InfoFileOld DataFileOld

[~, HName] = system('hostname');
HName = strsplit(HName);
if strcmp(HName{1}, 'will-linux')
    DataDir = '/home/will/Documents/data/OpTrap/';
elseif strcmp(HName{1},'bowtie')
    DataDir = '/home/ppxwh2/Documents/data/OpTrap/';
end

InfosDir = [DataDir 'infos/'];

if nargin >= 4
    ImstackOnly = varargin{1};
else 
    ImstackOnly = false;
end
if nargin >= 5 && ~isempty(varargin{2})
    if strncmp('_',varargin{2},1)
        InfoPrefix = varargin{2};
    else
        InfoPrefix = ['_' varargin{2}];
    end
else 
    InfoPrefix = '';
    warning('Default filenames were changed - if you want info_seg... use InfoPrefix')
end
if nargin >= 6 && ~isempty(varargin{3})
    if strncmp('_',varargin{3},1)
        InfoSuffix = varargin{3};
    else
        InfoSuffix = ['_' varargin{3}];
    end
else 
    InfoSuffix = '';
end
if nargin >= 7
    SimpleMatName = varargin{4};
else
    SimpleMatName = false;
end

% 1) Get filenames for data and info
if strcmp(CellType,'hela')
    InfoFile = {[InfosDir 'info-hela_ctrl_s_020_tr_70_' Num InfoSuffix '.mat']};
    DataFile = {'0610/Deformation','ctrl','020','70',0,Num};
    Loader = 'load_imstack';
elseif strcmp(CellType, 'helaS3')
    
    infoFile = {[InfosDir 'info-hela_' Set Num InfoSuffix '.mat']};
    DataFile = {'2020_11_19',
elseif strcmp(Set,'Jenna')
    InfoFile = {strcat(InfosDir,'info_Jenna_test_no_erode',Num,InfoSuffix,'.mat')};
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
    if ~SimpleMatName
        InfoFile = {[InfosDir 'info_reduced' InfoPrefix '_' CellType '_' Set matName InfoSuffix '.mat']};
    else
        InfoFile = {[InfosDir 'info_reduced' InfoPrefix '_' CellType '_' Set '_' Num InfoSuffix '.mat']};
    end
    DataFile = {[DataDir '2017_10_movies-from-aishah/'...
        CellType '/' CellType '_' Set '/' matName(2:end) '.avi']};
    Loader = 'avi_to_imstack';
elseif strcmp(CellType,'LS174T')
    switch Set
        case 'hypoxia'
            Date = '210717';
        case 'normoxia'
            Date = '200717';
    end
    Loader = 'avi_to_imstack';
    DataFile = {[DataDir '2017_10_movies-from-aishah/LS174T/' Date '_' Num '_LS174T_' Set '_1.avi']};
    InfoFile = {[InfosDir 'info_reduced' InfoPrefix '_LS174T_' Set '_' Date '_' Num '_LS174T_' Set '_1' InfoSuffix '.mat']};
elseif strcmp(CellType,'MV411')
    switch Set
        case 'normoxia'
            if str2double(Num) <= 10
                FName = ['100717_' Num '_ ' CellType '_1.avi'];
            elseif ~strcmp(Num, '25')
                FName = [CellType '_' Num '_0.020mms-1_1.avi'];
            else
                FName = [CellType '_' Num '_0.020mms-1_2.avi'];
            end
        case 'with_drugs'
            FName = ['180717_' Num '_' CellType '_0.020mms-1_1.avi'];
            if strcmp(Num,'12')
                FName(end-4) = '2';
            end
    end
    InfoFile = {[InfosDir 'info_reduced' InfoPrefix '_' CellType '_' Set '_' FName(1:end-4) InfoSuffix '.mat']};
    DataFile = {[DataDir '2017_10_movies-from-aishah/MV411/MV411_' Set '/' FName]};
    Loader = 'avi_to_imstack';
else
    error(['Cell type ' CellType ' not recognised'])
end

% 2) Load it if the last file was different to this
if ~strcmp(InfoFileOld, InfoFile) || ~strcmp(DataFileOld, DataFile) || ImstackOnly
    Imstack = eval([Loader '(DataFile{:});']);
    if ~ImstackOnly
        S = load(InfoFile{:});
        % 3) Extract info and meta for output
        info = S.info;
        meta = S.meta;
        % 4) Update trackers
        InfoFileOld = InfoFile;
    else
        info = [];
        meta = [];
    end
    DataFileOld = DataFile;
end

if nargout >= 1
    varargout{1} = Imstack;
end
if nargout >= 2
    varargout{2} = info;
end
if nargout >= 3
    varargout{3} = meta;
end

end