%% Show unwrapped result with errorbars
% An attempt at characterising errors from unwrap analysis compared to the
% "gold standard" errors from fitting to simulated data

%% Check data is correctly loaded
force_run_unwrap = false;
% Check correct info, meta and Imstack are loaded
compare_info_meta_imstack(info, meta, Imstack)

if isempty(whos('unwrapped')) || isempty(whos('FitErrs')) || isempty(whos('Ia')) || force_run_unwrap
    [u_fits, unwrapped, Ia, FitEqn, offset, FitErrs] = ...
        unwrap_cell_v2(Imstack, [info.mCentres] , repmat(100,1,size(Imstack{1},1)),'sc_up',1.8,'ifNaN','mean','sc_down',0.35);
end
N_frames = size(Imstack{1},1);

tol = 0.15;
idxa = Ia > (1-tol)*median(Ia,2) & Ia < (1+tol)*median(Ia,2); % Indexes for values included in fitting

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
errorbar([info.uTaylorParameter],DErrs,'.')
%% Compare these errors to simulated errors
% Take the dataset, bin it into 10 bands which will all have the same
% relative error. Calculate relative errors by simulating data for each
% deformation band, fitting simulated data and using error against ground
% truth as assumed error on actual result.

% need to rejig this so edges is passed into SimulateUnwrapFitting
%[N, Edges] = histcounts([info.uTaylorParameter],Dvals);
[Errors, Dvals] = SimulateUnwrapFitting(Ia);

Banding = sum([info.uTaylorParameter] > Dvals',1);

figure(88)
clf
hold on
errorbar([info.uTaylorParameter],Errors(Banding).*[info.uTaylorParameter])

%% Calculate error from residuals
% This produces errors compararable to standard deviation method - which
% makes sense if they have the same source
Thetas = repmat((1:360)',1,N_frames);

Residuals = squeeze(Ia) - FitEqn(u_fits(1,:), u_fits(2,:), u_fits(3,:),Thetas);
RMS = sqrt(mean(Residuals.^2,1));