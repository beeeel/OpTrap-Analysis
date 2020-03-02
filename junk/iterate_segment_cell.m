%% Optimise cell segmentation
% Use fminsearch to find the optimal settings for segment_cell_v2
% I don't know if this will work
clear all
%% Load datasets
dir = '0514';
treat = 'cytd';
speed = '014';
trap = '70';
spots = '3';
Imstack = load_imstack(dir, treat, speed, trap, spots, 0);
%% Use fminsearch for optimisation
opts = optimset('Display','iter','MaxFunEvals', 1000);
val_fun = @(par) segment_cell_v2_opt(Imstack, 'Area', par(1), par(2), par(3));
    % Options are: sigma, Ethresh, Fthresh
OptPar = fminsearch(val_fun, [1, 0.5473, 0.5739], opts); 
% Started with 2, 0.5, 51, 0.4: Opt stopped at 1, 0.5473, 53, 0.5739
% Modified segment_cell_v2_opt to use windowsize 53 fixed.
% Started with 1, 0.5473, 0.5739: Opt stopped at 1.0149, 0.5456, 0.5551