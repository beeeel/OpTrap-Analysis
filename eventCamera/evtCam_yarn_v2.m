%% Settings
evtCamDir = '~/Documents/data/EventCam/';
subDir = 'Tanias_data/';

sigmaFactor = 250;
nTimes = 2e5;

% idxPreExp = '1:numel(cdEvents.x)';
% idxPostExp = 'round(dT:dT:(length(cdEvents.ts)-dT))';
%% Prep 

%% Load data
addpath('~/Documents/Analysis/prophesee-matlab-scripts');   % Event cam loading scripts

dirList = dir([evtCamDir subDir]);
ls([evtCamDir subDir])

allTime = tic;
count = 1;
% fileBase = 'TrappedBead_3'; % 1, 2, or 3
for fIdx = [4:5 9]%:length(dirList)
    if length(dirList(fIdx).name) > 4
        if strcmp(dirList(fIdx).name(end-3:end), '.dat')
            fName = dirList(fIdx).name
            oneTime = tic;

            cdEvents = load_cd_events([evtCamDir subDir fName])
            
            %% Regularize time data for reasons
            sigma = sigmaFactor*range(cdEvents.ts)./nTimes;
            gaussWdth = 2*ceil(2*sigma)+1;
            
            nEvents = numel(cdEvents.ts);
            dT = floor(range(cdEvents.ts)/nTimes);
            evtRate = nEvents/range(cdEvents.ts);
            
            fprintf('Got %.2E events\t\t Range %is\t Mean rate %iHz\n',nEvents, range(cdEvents.ts)*1e-6, nEvents/range(cdEvents.ts)*1e-6)
            fprintf('Reconstruction dT %ius\t Mean events per point %.f\n', dT, gaussWdth*evtRate)
            
            times = round(0:dT:(max(cdEvents.ts)-dT-1))+min(cdEvents.ts);
            times2 = round(0:dT:(max(cdEvents.ts)-dT-1))+min(cdEvents.ts)+gaussWdth;
            
            tIdxs = zeros(length(times),1);
            tIdxs2 = zeros(length(times),1);
            
            residuals = zeros(length(times),1);
            residuals2 = zeros(length(times),1);
            
            %             tIdxs = tI{1};
%             tIdxs2 = tI{2};
%             idx = 0;
%             idxs = 1:2e5;
%             tic
%             while idx < length(times)
%                 idx = idx + 1;
%                 
%                 % You can go your own waaaayyy
% %                 [residuals(idx), tIdxs(idx)] = min(abs(cdEvents.ts(idxs) - times(idx)));
%                 [residuals(idx), tIdxs(idx)] = min(abs(cdEvents.ts - times(idx)));
% 
%                 
%                 % Needed better check for idxs going beyond max
% %                 while idxs(end) < nEvents  && residuals(idx) == abs(cdEvents.ts(idxs(end)) - times(idx))
% %                     idxs = idxs + min(15e3,nEvents-idxs(end));
% %                     [residuals(idx), tIdxs(idx)] = min(abs(cdEvents.ts(idxs) - times(idx)));
% %                 end
%                 tIdxs(idx) = tIdxs(idx) + idxs(1) - 1;
%             end
%             toc
%             disp('timed')
                
            % Get indexes for raw data
            tI = cell(2,1);
            ts = cdEvents.ts;
            alltimes = [times' times2'];
            tic
            for t = 1:2
                time = alltimes(:,t);
                tIdxs = NaN(length(time),1);
                idx = 0;
                while idx < length(time)
                    idx = idx + 1;
                    tmp = ts - time(idx);
                    tmpIdx = sum(tmp <= 0);
                    tIdxs(idx) = tmpIdx;
                    residuals = tmp(tmpIdx);
                end
                tI{t} = tIdxs;
                toc
            end
            idx = 0;
            idxs = 1:2e5;
            while idx < length(times2)
                idx = idx + 1;
                
                % You can go your own waaaayyy
%                 [residuals2(idx), tIdxs2(idx)] = min(abs(cdEvents.ts(idxs) - times2(idx)));
                [residuals2(idx), tIdxs2(idx)] = min(abs(cdEvents.ts - times2(idx)));
                
%                 while idxs(end) < nEvents && residuals2(idx) == abs(cdEvents.ts(idxs(end)) - times2(idx))
%                     idxs = idxs + min(15e3,nEvents-idxs(end));
%                     [residuals2(idx), tIdxs2(idx)] = min(abs(cdEvents.ts(idxs) - times2(idx)));
%                 end
                tIdxs2(idx) = tIdxs2(idx) + idxs(1) - 1;
            end
            disp('time2''d')
            
            eventsGrouped = cell(length(tIdxs),1);
            for idx = 1:length(tIdxs)
                idxs = tIdxs(idx):tIdxs2(idx);
                %     eventsGrouped{idx} = [cdEvents.ts(idxs) cdEvents.x(idxs) cdEvents.y(idxs) cdEvents.p(idxs)];
                eventsGrouped{idx} = [cdEvents.ts(idxs) cdEvents.x(idxs) cdEvents.y(idxs)];
            end
            disp('grouped')
            
            % Apply filter in parallel
            xGauss = zeros(length(tIdxs),1);
            yGauss = zeros(length(tIdxs),1);
            tGauss = zeros(length(tIdxs),1);
            nGauss = zeros(length(tIdxs),1);
            
            parfor idx = 1:length(tIdxs)
                data = eventsGrouped{idx};
                
                mu = (times2(idx) + times(idx))*0.5;
                gauss = (1/sqrt(2*pi*sigma.^2)) * exp(- (data(:,1) - mu).^2/sigma.^2);
                normFactor = sum(gauss);
                
                tGauss(idx) = sum(data(:,1).*gauss)./normFactor;
                xGauss(idx) = sum(data(:,2).*gauss)./normFactor;
                yGauss(idx) = sum(data(:,3).*gauss)./normFactor;
                nGauss = size(data,1);
            end
            disp('done!')
            
            %% MSDs
            tracks = {[1e-6*tGauss, xGauss], [1e-6*tGauss, yGauss]};
            
            msd = msdanalyzer(1, 'px', 's','log');
            msd = msd.addAll(tracks);
            msd = msd.computeMSD;
            
            figure(1)
            clf
            ax = gca;
            msd.plotMSD;
            ax.XAxis.Scale = 'log';
            ax.YAxis.Scale = 'log';
            legend('X','Y','Location','best')
            
            %% Save
            save([evtCamDir subDir 'MSDs/' fName(1:end-4) '_msd.mat'], 'msd', 'fName','xGauss','yGauss','tGauss','nGauss')
            fprintf('Completed %i in %0.1f\t\t at %0.1f\n',count, toc(oneTime), toc(allTime));
            count = count + 1;
        end
    end
end

cmd = ['tar -C ' evtCamDir subDir ' -cvf MSDs/msds.tar *.mat'];
system(cmd)
