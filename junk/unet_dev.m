imageSize = [repmat(2^ceil(log2(max(size(Imstack{1}{1,1})))),1,2), 1];
numClasses = 2;
lgraph = unetLayers(imageSize,numClasses);

imds = imageDatastore('~/png/ML/Samples','IncludeSubFolders',true,'ReadSize',1e3);
% Need to create pixel label datastore to do segmentation, or having
% training data some other way
numTrainingFiles = 0.75 * length(imds.Files);

[imdsTrain,imdsTest] = splitEachLabel(imds,numTrainingFiles,'randomize'); % Images need labels for this

options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.001, ...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(imdsTrain, lgraph, options);