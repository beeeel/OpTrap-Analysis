%% Plot data from bead_processing_accumulator
%
open ./bead_processing_accumulator_v0.m
%% First load it
clear 
close all
run ./bead_processing_accumulator_v0.m
fprintf('Loaded days: \t')
disp(dayDirs)

daysLegCell = cell(1,length(dayDirs));
for dIdx = 1:length(dayDirs)
    tmp = strsplit(dayDirs{dIdx},{'_' '/'});
    daysLegCell{dIdx} = strjoin(tmp(3:-1:1), '/');
    
    fprintf('\n\n  Date: %s \n  Variable ''accumulated'' contains: \n', daysLegCell{dIdx})
    disp(accumulated{dIdx})
end
%% Plot aesthetics
% Choose which field from accumulated to plot
fieldName = 'stiffraw';
% Set colours and marks
lineStyles = {'--' ':' '-.' '-'};
colours = {'b' 'k' 'm' [0 0.7 0]};
marks = '^oxv';
markSz = 5;
%% When were the measurements taken?
figure(200)
clf
hold on

% First get the colours for the legend
cCount = 1;
for dIdx = 1:length(accumulated)
    plot(thymes{dIdx}{1}(1), cCount, marks(dIdx), 'Color', colours{dIdx}, 'MarkerSize', markSz)
    cCount = cCount + size(accumulated{dIdx},2);
end

% Then do all the data
cCount = 1;
for dIdx = 1:length(accumulated)
    for cIdx = 1:size(accumulated{dIdx},2)
        plot(thymes{dIdx}{cIdx}, cCount * ones(1,length(thymes{dIdx}{cIdx})), marks(dIdx), 'Color', colours{dIdx}, 'MarkerSize', markSz)
        cCount = cCount + 1;
    end
end

legend(daysLegCell)
% legend('Latrunculin added','Location','best')
% legend('10% Water added','Location','best')

ylim([0.5 cCount+0.5])
xlabel('Time after seeding (min)')
ylabel('Cell number')
title('Distribution of measurements in time')
%% Low-frequency elastic measurements
figure(201)
clf

XY = 'XY';
ab = {'(a)','(b)'};

% Do for X and Y
for plt = 1:2
    ax = subplot(1,2,plt);
    ax.YScale = 'log';
    hold on
    
    % First get the colours for the legend
    for dIdx = 1:length(accumulated)
        tmp = accumulated{dIdx}{1,1}(1);
        % Convert N/m to μN/m
        stiff = 1e6 * tmp.(fieldName);
        plot(thymes{dIdx}{1}(1), stiff(plt), 'LineStyle', lineStyles{dIdx}, ...
            'Color', colours{dIdx}, 'MarkerSize', markSz, 'Marker', marks(dIdx))
    end
    
    % Then do all the data
    for dIdx = 1%:length(accumulated)
        for cIdx = 2%1:size(accumulated{dIdx},2)
            tmp = accumulated{dIdx}{1,cIdx};
            % Here stiff goes [X1 Y1 X2 Y2 ...]
            stiff = 1e6 * [tmp.(fieldName)];
            plot(thymes{dIdx}{cIdx}, stiff(plt:2:end), 'LineStyle', lineStyles{dIdx}, ...
                'Color', colours{dIdx}, 'MarkerSize', markSz, 'Marker', marks(dIdx))
        end
    end
    
    legend(daysLegCell,'Location','best')
    title(sprintf('Low frequency elastics \n%s direction', XY(plt)))
%     title(sprintf('%s',ab{plt}))
    xlabel('Time after seeding (min)')
    ylabel('Stiffness βG''_0 (μN/m)')
end
%% XY Ratio
figure(202)
clf
hold on

% Assume to be using the same plot aesthetics as above

% First get the colours for the legend
for dIdx = 1:length(accumulated)
    tmp = accumulated{dIdx}{1,1}(1);
    stiff = tmp.(fieldName)(1) / tmp.(fieldName)(2);
    plot(thymes{dIdx}{1}(1), stiff, 'LineStyle', lineStyles{dIdx}, ...
        'Color', colours{dIdx}, 'MarkerSize', markSz, 'Marker', marks(dIdx))
end

% Then do all the data
for dIdx = 1:length(accumulated)
    for cIdx = 1:size(accumulated{dIdx},2)
        tmp = accumulated{dIdx}{1,cIdx};
        % Here stiff goes [X1 Y1 X2 Y2 ...]
        stiff = [tmp.(fieldName)];
        stiff = stiff(1:2:end) ./ stiff(2:2:end);
        plot(thymes{dIdx}{cIdx}, stiff, 'LineStyle', lineStyles{dIdx}, ...
            'Color', colours{dIdx}, 'MarkerSize', markSz, 'Marker', marks(dIdx))
    end
end

ax = gca;
ax.YScale = 'log';

plot(xlim, [1 1], '-', 'Color', [1 1 1] .* 0.6, 'LineWidth', 2)

legend(daysLegCell,'Location','best')
title('Low frequency elastics ratio X ÷ Y')
xlabel('Time after seeding (min)')
ylabel('Stiffness βG''_0 ratio ( X ÷ Y )')
%% MSDs
dIdx = 1;
cIdx = 1;

% colours = {[0.3 0.5 0.9] [0.2 0.4 0.8] [0.1 0.3 0.7] [};
tits = {'Radial','Tangential'};

figure(203)
clf

msdObj = [accumulated{dIdx}{1,cIdx}.msd];
ts = accumulated{dIdx}{2,cIdx};

colormap cool
colourmap = colormap;
colour = colourmap(ceil(size(colourmap,1)*(ts-ts(1)+1e-6)/(ts(end) - ts(1)+1e-6)),:);

for plt = 1:2
    ax = subplot(1,2,plt);
    hold on
    for rep = 1:length(msdObj)
        plot(msdObj(rep).msd{plt}(:,1), msdObj(rep).msd{plt}(:,2), ...
            'Color',colour(rep,:))
    end
    ax.XScale = 'log';
    ax.YScale = 'log';
    xlabel('Delay time \tau (s)')
    ylabel('Mean-squared displacement (\mu m^2)')
    title(sprintf('%s direction',tits{plt}))
    axis equal
end
colorbar 
ax.CLim = [ts(1) ts(end)];
%% NMSDs
dIdx = 1;
cIdx = 2;

sf = [22 15];

fSz = 16;
ab = {'(a)','(b)'};

figure(204)
clf

msdObj = [accumulated{dIdx}{1,cIdx}.msd];
ts = accumulated{dIdx}{2,cIdx};

colormap cool
colourmap = colormap;
colour = colourmap(ceil(size(colourmap,1)*(ts-ts(1)+1e-6)/(ts(end) - ts(1)+1e-6)),:);

for plt = 1:2
    ax = subplot(1,2,plt);
    hold on
    X = logspace(-3.1, -1, 10);
    plot(X, sf(1).*X,'k--','LineWidth',2)
    plot(X, sf(2).*X.^0.75,'r--','LineWidth',2)
    for rep = 1:length(msdObj)
        normFac = var(msdObj(rep).tracks{plt}(:,2));
        plot(msdObj(rep).msd{plt}(:,1), msdObj(rep).msd{plt}(:,2)./normFac, ...
            'Color',colour(rep,:), 'LineWidth', 2)
    end
    ax.XScale = 'log';
    ax.YScale = 'log';
    ax.FontSize = fSz;
    xlabel('Delay time \tau (s)')
    ylabel('\Pi (\tau)')
    title(ab{plt})
    legend('\alpha\tau^1','\alpha\tau^{0.75}','Location','southeast')
end
colorbar 
ax.CLim = [ts(1) ts(end)];
%% NPAF

dIdx = 2;
cIdx = 3;

msd = [accumulated{dIdx}{1,cIdx}.msd];

for idx = 1:4
    nmsdX = msd(idx).msd{1}(:,2) ./ var(msd(idx).tracks{1}(:,2));
    NPAFX{idx} = 1 - nmsdX;
    
    nmsdY = msd(idx).msd{2}(:,2) ./ var(msd(idx).tracks{2}(:,2));
    NPAFY{idx} = 1 - nmsdY;
end
%% plot it
figure(205)
clf

for idx = 1:4
    ax = subplot(1,2,1);
    plot(msd(idx).msd{1}(:,1), NPAFX{idx})
    ax.XScale = 'log';
    ax = subplot(1,2,2);
    plot(msd(idx).msd{2}(:,1), NPAFY{idx})
    ax.XScale = 'log';
end
%% FFT it
figure(206)
clf
hold on
for idx = 1:4
    ax = subplot(1,2,1);
    [w, X] = fft_scaled( msd(idx).msd{1}(:,1), NPAFX{idx}, true, ax);
    
    ax = subplot(1,2,2);
    [w, X] = fft_scaled( msd(idx).msd{2}(:,1), NPAFY{idx}, true, ax);
end
%% Average first MSD for each cell
allTracksX = {};
allTracksY = {};
idx = 1;
nT = 5e4;
for dIdx = 1:length(accumulated)
    for cIdx = 1:length(accumulated{dIdx})
        for tIdx = 0:(length(accumulated{dIdx}{1,cIdx}(1).msd.tracks{1})/nT-1)
            allTracksX{idx} = accumulated{dIdx}{1,cIdx}(1).msd.tracks{1}(tIdx*nT+1:(1+tIdx)*nT,:);
            allTracksY{idx} = accumulated{dIdx}{1,cIdx}(1).msd.tracks{2}(tIdx*nT+1:(1+tIdx)*nT,:);
            idx = idx + 1;
        end
    end
end

%% 
msdX = msdanalyzer(1, 'um', 's');
msdY = msdanalyzer(1, 'um', 's');
msdX = msdX.addAll(allTracksX);
msdY = msdY.addAll(allTracksY);
msdX = msdX.computeMSD;
msdY = msdY.computeMSD;
%%
% 2021_07_06 : Calculated msdX took 30 minutes. Saved as "MeanMSD_X.mat".
% Also saved MeanMSD_Y.mat
% The .mat files ↑ probably don't have the data I wanted. There's a .fig
% which will contain it as XData and YData for 2 axis, mean_MSD_firstDataEachCell_50k.fig
figure(1)
ax = subplot(2,1,1);
msdX.plotMeanMSD;
ax.XScale = 'log';
ax.YScale = 'log';

ax = subplot(2,1,2);
msdY.plotMeanMSD;
ax.XScale = 'log';
ax.YScale = 'log';
%%
fSz = 16;
ax = subplot(2,1,1);
title('X direction');
ax.FontSize = fSz;
ax = subplot(2,1,2);
title('Y direction');
ax.FontSize = fSz;
%%
%% MSDs with fits
dIdx = 1;
cIdx = 1;

% tRange = [2 70];
tRange = [0.004 0.15];
% colours = {[0.3 0.5 0.9] [0.2 0.4 0.8] [0.1 0.3 0.7] [};
tits = {'Radial','Tangential'};

figure(203)
clf

msdObj = [accumulated{dIdx}{1,cIdx}.msd];
ts = accumulated{dIdx}{2,cIdx};

colormap cool
colourmap = colormap;
colour = colourmap(ceil(size(colourmap,1)*(ts-ts(1)+1e-6)/(ts(end) - ts(1)+1e-6)),:);

for plt = 1:2
    ax = subplot(1,2,plt);
    hold on
    ah = annotation('textbox', [(0.35+(plt-1)*0.15) 0.85 0.15 0.15]);
    ah.BackgroundColor = 'white';
    for rep = 1:length(msdObj)
        plot(msdObj(rep).msd{plt}(:,1), msdObj(rep).msd{plt}(:,2), ...
            'Color',colour(rep,:))
        msdIdx = find(msdObj(rep).msd{plt}(:,1) > tRange(1), 1) ...
            : ( find(msdObj(rep).msd{plt}(:,1) > tRange(2), 1) - 1 );
        tauData = msdObj(rep).msd{plt}(msdIdx,1);
        msdData = msdObj(rep).msd{plt}(msdIdx,2);
        fo = fit(log(tauData), log(msdData), 'poly1');
        plot(tauData, exp(fo.p1 * log(tauData))*exp(fo.p2), 'k:', 'LineWidth', 1.5)
        ah.String{rep} = sprintf('t = %i, MSD \\propto τ^{%0.2g}', ts(rep), fo.p1);
    end
    ax.XScale = 'log';
    ax.YScale = 'log';
    xlabel('Delay time \tau (s)')
    ylabel('Mean-squared displacement (\mu m^2)')
    title(sprintf('%s direction',tits{plt}))
    axis equal
end
colorbar 
ax.CLim = [ts(1) ts(end)];