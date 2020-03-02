%% Tune cell finding parameters to best find cells
% Use fminsearch on a modified version of find_cell to optimise 
% Highly sensitive to initial conditions, doesn't always find sensible
% results
%% Load datasets
dir = '0610/Deformation';
treat = 'ctrl';
speed = '020';
trap = '70';
spots = 0;
rep = '1';
Imstack = load_imstack(dir, treat, speed, trap, spots, rep);
%% Use fminsearch for optimisation
opts = optimset('Display','iter','MaxFunEval', 400);
val_fun = @(par) find_cell_v2_opt(Imstack, par(1), par(2), par(3), par(4));
start = [10, 250, 0.9, 2];
lb = [1, 0, 0, 1];
ub = [50, 1000, 1, 5];
[OptPar, fval, exitflag, output] = patternsearch(val_fun, start, [],[],[],[], lb, ub);
% Started at 10, 250, 0.9, 2: ended near 10, 10. 

% results on cytd sets:
% Started at 20, 20: ended near 20, 21 - gave shit results
% Started at 15, 15: ended near 20, 16