function data = bead_loadData(data)
%% data = loadBeadData(data)
% Load centres, times, images and metadata for a given dataset

data.raw.xCentresPx = byteStreamToDouble([data.dirPath '/XCentres.dat']);
data.raw.yCentresPx = byteStreamToDouble([data.dirPath '/YCentres.dat']);
data.raw.timeVecMs  = byteStreamToDouble([data.dirPath '/Times.dat']);
data.Imstack  = bfopen([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
metadata = fileread([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default_metadata.txt']);
data.metadata = jsondecode(metadata);
data.nPoints = length(data.raw.xCentresPx);