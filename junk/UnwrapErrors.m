%% Show unwrapped result with errorbars
% An attempt at characterising errors from unwrap analysis compared to the
% "gold standard" errors from fitting to simulated data

%% Check data is correctly loaded
force_run_unwrap = false;
% Check correct info, meta and Imstack are loaded

if ~isempty(whos('info')) && ~isempty(whos('meta')) && ~isempty(whos('Imstack'))
    compare_info_meta_imstack(info, meta, Imstack)
end

if force_run_unwrap || isempty(whos('Ia')) || isempty(whos('unwrapped'))
    [u_fits, ~, Ia, FitEqn, offset, FitErrs] = ...
        unwrap_cell_v2(Imstack, [info.mCentres] , repmat(100,1,size(Imstack{1},1)),'sc_up',1.8,'ifNaN','mean','sc_down',0.35); %#ok<UNRCH>
end
N_frames = size(Imstack{1},1);

tol = 0.15;
idxa = Ia > (1-tol)*median(Ia,2) & Ia < (1+tol)*median(Ia,2); % Indexes for values included in fitting

%% Plotting variables
FSize = 16;

Tdata = linspace(0,9.99,size(Imstack{1},1));

%% Calculate error from standard deviation of relaxed cell edge position
% This doesn't produce big enough errors 
Sum = sum(u_fits(1:2,:),1);
Diff = u_fits(1,:) - u_fits(2,:);
StDev = reshape(std(Ia,0,2),[],N_frames);

DRelErrs = ((StDev./u_fits(1,:)).*((1./Sum) - Diff./Sum.^2)).^2 + ...
    ((StDev./u_fits(2,:)).*((-1./Sum) - Diff./Sum.^2)).^2;
DErrs = DRelErrs.* [info.uTaylorParameter];

figure(87)
clf
hold on
errorbar(Tdata, [info.uTaylorParameter],DErrs,'.')
xlabel('Time (s)','FontSize',FSize)
ylabel('Deformation','FontSize',FSize)
title({'Deformation with effors from standard' 'deviation of relaxed cell edge position'},...
    'FontSize',FSize)

%% Compare these errors to simulated errors
% Take the dataset, bin it into 10 bands which will all have the same
% relative error. Calculate relative errors by simulating data for each
% deformation band, fitting simulated data and using error against ground
% truth as assumed error on actual result.

NBins = 100;
[N, Edges] = histcounts([info.uTaylorParameter],NBins);
if Edges(1) == 0; Edges = Edges(2:end); end
[Errors, Dvals] = SimulateUnwrapFitting(Ia,Edges);

Banding = sum([info.uTaylorParameter] > Dvals',1);
Banding(Banding==0) = 1;

%%

figure(88)
clf
hold on
Line = errorbar(Tdata,[info.uTaylorParameter],Errors(Banding).*[info.uTaylorParameter]);
Line.Marker = 'x';
Line.MarkerEdgeColor = 'r';
title({'Deformation with errors from simulated data' 'using noise similar to first 100 frames'},'FontSize',FSize)
xlabel('Time (s)','FontSize',FSize)
ylabel('Deformation','FontSize',FSize)
%% Calculate error from residuals
% This produces errors compararable to standard deviation method - which
% makes sense if they have the same source
Thetas = repmat((1:360)',1,N_frames);

Residuals = squeeze(Ia) - FitEqn(u_fits(1,:), u_fits(2,:), u_fits(3,:),Thetas);
RMS = sqrt(mean(Residuals.^2,1));

%% Augment data with rotations
Frames = 1:10:1000;
NRotations = 10;

Rotations = linspace(0, 2*pi, NRotations);
RCentres = zeros(2,length(Frames)*length(Rotations));
FrCount = 1;
for frame = Frames
    for rot = Rotations
        AugStack{1}{FrCount,1} = imrotate(Imstack{1}{frame,1},180*rot/pi,'bilinear','crop');
        AugStack{1}{FrCount,2} = [Imstack{1}{frame,2} ' rotated ' num2str(rot)];
        RCentres(:,FrCount) = info(frame).mCentres;
        FrCount = FrCount + 1;
    end
end

[aug_fits, ~, aug_Ia, ~, ~, FitErrs] = ...
        unwrap_cell_v2(AugStack, RCentres , repmat(100,1,size(AugStack{1},1)),'sc_up',1.8,'ifNaN','mean','sc_down',0.35);
%%
XData = repmat(Rotations,1,length(Frames)) + reshape(2*pi*repmat(0:length(Frames)-1,length(Rotations),1),1,[]);
figure(89)
clf
hold on
Line = errorbar(XData,aug_fits(3,:),FitErrs(3,:));
plot(2*pi*reshape([1; 1].* (1:length(Frames)),1,[]),0.5*pi*reshape([-1; 1; 1; -1].* ones(1,length(Frames)/2),1,[]),'r--')

figure(90)
clf
hold on
plot(XData,aug_fits(1,:))
plot(XData,aug_fits(2,:))
%% Compare fitting with and without rotations

AugD = Fits2Ds(aug_fits);
FitD = Fits2Ds(u_fits);
AugErrs = AugD-repmat(FitD(Frames),1,NRotations);
PltErrs = mean(reshape(AugErrs,1,[],NRotations),3);

%AugErrs = aug_fits - repmat(u_fits(:,Frames),1,NRotations);
%PltErrs = squeeze(mean(abs(reshape(AugErrs,3,[],NRotations)),3));

%PltFits = squeeze(mean(reshape(aug_fits,3,[],NRotations),3));
%%
figure(91)
clf
hold on
errorbar(Frames,FitD(Frames),PltErrs(1,:))
%%
function Ds = Fits2Ds(fits)
Ds = (fits(1,:) - fits(2,:))./(fits(1,:) + fits(2,:));
end

function R = Rot2D(angle)
    R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
end