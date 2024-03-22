function data = bead_loadData(data, varargin)
%% data = loadBeadData(data/dirPath, [loadImages, skipSuffixes])
% Load centres, times, images and metadata for a given dataset

% Make the struct
if ~isstruct(data)
    if strcmp(data(1), '/')
        dP = data;
        fN = strsplit(data, '/');
        fN = fN{end};
    else
        dP = [pwd '/' data];
        fN = data;
    end
    data = struct('fName', fN, 'dirPath',dP);
end

% Load opts
data.opts = bead_loadOpts(data);

% Handle inputs
loadImages = true;
if nargin > 1
    loadImages = varargin{1};
end
skipSuffixes = {};
if nargin > 2
    skipSuffixes = varargin{2};
end

if ~isfield(data.opts, 'skipSuffixes')
    data.opts.skipSuffixes = skipSuffixes;
else
    if ~isempty(skipSuffixes)
        warning('skipSuffixes overridden by opts file')
    end
    warning('I removed a line of code because it seemed stupid, now that line didn''t load')
%     data.opts.skipSuffixes = strsplit(data.opts.skipSuffixes,',');
end

% New method: Find anything with a suffix and then record it and what
% suffix
if ~isfield(data, 'raw') ||  ~isfield(data.raw, 'suffixes')
    data.raw.suffixes = {};
end
errMsg = '';

if ~exist(data.dirPath,'dir')
    error('Folder %s does not exist!',data.dirPath)
end

dl = dir(data.dirPath);
il = dl(endsWith({dl.name}, '.dat') & startsWith({dl.name}, {'I'}));
dl = dl(endsWith({dl.name}, '.dat') & startsWith({dl.name}, {'X', 'Y','T','Times'}));
nx = sum(startsWith({dl.name}, 'X'));
ny = sum(startsWith({dl.name}, 'Y'));
nt = sum(startsWith({dl.name}, {'T','Times'}));

checkSizes;

tl = dl(startsWith({dl.name}, {'T','Times'}));
[~, tidx] = min([tl.bytes]);
data.raw.timeVecMs  = byteStreamToDouble([data.dirPath '/' tl(tidx).name]);
nP = length(data.raw.timeVecMs);

nNans = zeros(nx,1);
xdx = 1;
ydx = 1;
for idx = 1:length(dl)
    fName = sprintf('%s/%s', data.dirPath, dl(idx).name);
    suff = strsplit(dl(idx).name, {'X', 'Y', '.dat'});
    suff = suff{end-1};
    if any(strcmp(data.opts.skipSuffixes, suff))
%         warning('Skipping file %s because skipSuffixes',dl(idx).name)
    else
        % Screen for NaNs
        dP = byteStreamToDouble(fName);
        if length(dP) ~= nP
            warning('Skipping file %s with %i points (Times has %i points)', fName, length(dP), nP)
        else
            nNans(xdx) = warnNaN(dP, fName);
            % This next line might cause errors down the line but tbh I
            % should be doing it.
%             dP(isnan(dP)) = mean(dP(~isnan(dP)));
            if startsWith(dl(idx).name, 'X')
                data.raw.xCentresPx(xdx,:) = dP;
                if length(data.raw.suffixes)<xdx || strcmp(data.raw.suffixes{xdx}, suff)
                    data.raw.suffixes{xdx} = suff;
                else
                    error('Suffix mismatch like you said wouldn''t happen')
                end
                xdx = xdx + 1;
            elseif startsWith(dl(idx).name, 'Y')
                data.raw.yCentresPx(ydx,:) = dP;
                if length(data.raw.suffixes)<ydx || strcmp(data.raw.suffixes{ydx}, suff)
                    data.raw.suffixes{ydx} = suff;
                else
                    error('Suffix mismatch like you said wouldn''t happen')
                end
                ydx = ydx + 1;
            end
        end
    end
end
data.raw.nNans = nNans;

if ~isfield(data.opts, 'centresRow')
    if isfield(data.opts, 'thresh')
        data.opts.centresRow = find(endsWith(data.raw.suffixes, num2str(data.opts.thresh)),1);
        if isempty(data.opts.centresRow)
            error('Didn''t load data with thresh = %i... You knew you would want to do this eventually.', data.opts.thresh)
        end
    else
        data.opts.centresRow = 1:size(data.raw.xCentresPx,1);
    end
end

if exist([data.dirPath '/subWidth.dat'],'file')
    data.raw.subWidth = byteStreamToDouble([data.dirPath '/subWidth.dat']);
end

if isempty(il) 
    warning('No brightness (Z/I) data found')
else
    if isfield(data.opts, 'zthresh') 
        str = sprintf('I%%s_th%i.dat',data.opts.zthresh);
    elseif isfield(data.opts, 'thresh')
        str = sprintf('I%%s_th%i.dat',data.opts.thresh);
    else
        str = sprintf('I%%s_th0.dat');
    end
    for idx = 1:length(data.raw.suffixes)
        Ipath = [data.dirPath '/' sprintf(str,data.raw.suffixes{idx})];
        if exist(Ipath,'file')
            data.raw.dcAvg(idx,:) = byteStreamToDouble(Ipath);
        end
    end
end

if loadImages
    % More sensible method: If the folder exists, find all .tif files
    % within it. If there's only one, load it. Likewise for metadata.txt.
    if exist([data.dirPath '/images_and_metadata/'], 'dir')
        dirList = dir([data.dirPath '/images_and_metadata/']);
        imList = dirList(endsWith({dirList.name}, '.tif'));
        if isscalar(imList)
            fPath = strjoin({imList.folder imList.name}, '/');
            data.Imstack = bfopen(fPath);
        else
            error('Found multiple (or 0) TIFs in ROI folder')
        end
        txtList = dirList(endsWith({dirList.name}, 'metadata.txt'));
        
        if isscalar(txtList)
            fPath = strjoin({txtList.folder txtList.name}, '/');
            metadata = fileread(fPath);
            data.metadata = jsondecode(metadata);
            
            % Get ROI position
            roi = str2double( strsplit( data.metadata.FrameKey_0_0_0.ROI, '-' ) );
            data.opts.roi = roi;
        else
            error('Found multiple (or 0) files ending with metadata.txt in ROI folder')
        end
    else
        warning('Could not find ROI images or metadata')
    end
    
    if exist([data.dirPath '/full_images_and_metadata/'], 'dir')
        dirList = dir([data.dirPath '/full_images_and_metadata/']);
        dirList = dirList(endsWith({dirList.name}, 'tif'));
        if isscalar(dirList)
            fPath = strjoin({dirList.folder dirList.name}, '/');
            data.ImstackFullFoV  = bfopen(fPath);
        else
            error('Found multiple (or 0) TIFs in full FoV folder')
        end
    else
        warning('Could not find directory for full FoV images')
    end
else
    warning('Instructed to not load images')
end

data.nPoints = length(data.raw.xCentresPx);

function nNans = warnNaN(arr, name)
nNans = sum(isnan(arr));
if nNans
    warning(['Loaded and replaced ' num2str(nNans) ' from ' name])
end

end

function checkSizes

if nx ~= ny
    error('Found mismatched number of X (%i) and Y (%i) data', nx, ny)
elseif nt ~= 1 && mod(nx, nt) ~= 0
    error('Found mismatched number of T (%i) and XY (%i) data', nt, nx)
elseif all([nx ny nt] == 0)
    fprintf(errMsg);
    error('Could not find any centres data');
end

end

end
