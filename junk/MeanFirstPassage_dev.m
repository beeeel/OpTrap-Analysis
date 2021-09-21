%% Look at mean first passage time
% First passage time is the time taken to travel more than D distance from
% starting position. Average over all possible starting times.

% Set up a 2D logical operation
% column vector: distance of each observation from chosen starting position
% row vector: chosen distances, Ds.
% result: true if distance < Ds
% sum trues to find number of timesteps before passing distance D
% Place into running totals array
% Need sum, sum of squares, number of observations. 
% Calculate mean and std, convert into units of time.