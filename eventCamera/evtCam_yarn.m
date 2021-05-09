%% Settings
evtCamDir = '~/Documents/data/EventCam/';
subDir = 'Tanias_data/';

sigma = 500;
idxPreExp = '1:numel(cdEvents.x)';
idxPostExp = '1:10:numel(xGauss)';
%% Prep 
gaussSz = 2*ceil(2*sigma)+1;
gauss = 1:gaussSz;
gauss = exp(-(gauss-ceil(gaussSz/1)).^2./sigma.^2);
%% Load data
addpath('~/Documents/Analysis/prophesee-matlab-scripts');   % Event cam loading scripts
addpath('~/Documents/Analysis/tinevez-msdanalyzer-da0bdaa');% MSDanalyzer

dirList = dir([evtCamDir subDir]);
ls([evtCamDir subDir])

allTime = tic;
count = 1;
% fileBase = 'TrappedBead_3'; % 1, 2, or 3
for fIdx = 1:length(dirList)
    if length(dirList(fIdx).name) > 4
        if strcmp(dirList(fIdx).name(end-3:end), '.dat')
            fName = dirList(fIdx).name;
            oneTime = tic;

            cdEvents = load_cd_events([evtCamDir subDir fName])
            
            %% Gaussian filtering            
            idxPre = eval(idxPreExp);
            tmp = cdEvents.x(idxPre);
            xGauss = conv(tmp(cdEvents.p(idxPre) > 0), gauss, 'valid')./sum(gauss);
            tmp = cdEvents.y(idxPre);
            yGauss = conv(tmp(cdEvents.p(idxPre) > 0), gauss, 'valid')./sum(gauss);
            tmp = cdEvents.ts(idxPre);
            tGauss = conv(tmp(cdEvents.p(idxPre) > 0), gauss, 'valid')./sum(gauss);
            
            %% Subsampling to speed things up
            idxPost = eval(idxPostExp);
            xGauss = xGauss(idxPost);
            yGauss = yGauss(idxPost);
            tGauss = tGauss(idxPost);
            
            %% MSDs
            tracks = {[tGauss, xGauss], [tGauss, yGauss]};
            
            msd = msdanalyzer(1, 'px', 'us','original');
            msd = msd.addAll(tracks);
            msd = msd.computeMSDparfor;
            
            %% Save
            save([evtCamDir subDir fName(1:end-4) '_msd.mat'], 'msd', 'fName')
            fprintf('Completed %i in %0.1f\t\t at %0.1f\n',count, toc(oneTime), toc(allTime));
            count = count + 1;
        end
    end
end
