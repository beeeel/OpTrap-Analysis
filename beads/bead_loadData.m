function data = bead_loadData(data)
%% data = loadBeadData(data)
% Load centres, times, images and metadata for a given dataset

if exist([data.dirPath '/XCentres.dat'],'file')
    data.raw.xCentresPx = byteStreamToDouble([data.dirPath '/XCentres.dat']);
    data.raw.yCentresPx = byteStreamToDouble([data.dirPath '/YCentres.dat']);
elseif exist([data.dirPath '/X.dat'], 'file')
    data.raw.xCentresPx = byteStreamToDouble([data.dirPath '/X.dat']);
    data.raw.yCentresPx = byteStreamToDouble([data.dirPath '/Y.dat']);
elseif exist([data.dirPath '/Xl.dat'], 'file')
    data.raw.xCentresPx = byteStreamToDouble([data.dirPath '/Xl.dat']);
    data.raw.yCentresPx = byteStreamToDouble([data.dirPath '/Yl.dat']);
    data.raw.xCentresPx(2,:) = byteStreamToDouble([data.dirPath '/Xr.dat']);
    data.raw.yCentresPx(2,:) = byteStreamToDouble([data.dirPath '/Yr.dat']);
else
    disp(['File ' data.dirPath '/XCentres.dat does not exist'])
    disp(['File ' data.dirPath '/X.dat does not exist'])
    disp(['File ' data.dirPath '/Xl.dat does not exist'])
    error('Could not find any centres data');
end
if exist([data.dirPath '/I.dat'], 'file')
    data.raw.dcAvg = byteStreamToDouble([data.dirPath '/I.dat']);
end

if exist([data.dirPath '/Times.dat'], 'file')
    data.raw.timeVecMs  = byteStreamToDouble([data.dirPath '/Times.dat']);
elseif exist([data.dirPath '/T.dat'], 'file')
    data.raw.timeVecMs  = byteStreamToDouble([data.dirPath '/T.dat']);
else
    disp(['File ' data.dirPath '/Times.dat does not exist'])
    disp(['File ' data.dirPath '/T.dat does not exist'])
    error('Could not find any times data');
end

if exist([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif'], 'file')
    data.Imstack  = bfopen([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    metadata = fileread([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default_metadata.txt']);
    data.metadata = jsondecode(metadata);
else
    warning('Could not find ROI images or metadata')
end

if exist([data.dirPath '/full_images_and_metadata/full_images_and_metadata_MMStack_Default.ome.tif'], 'file')
    data.ImstackFullFoV  = bfopen([data.dirPath '/full_images_and_metadata/full_images_and_metadata_MMStack_Default.ome.tif']);
else
    warning('Could not find full FoV images')
end

data.nPoints = length(data.raw.xCentresPx);