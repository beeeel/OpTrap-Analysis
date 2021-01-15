function lT = calcLagTimes(nP, varargin)
%% Calculate lag times to use with allanvar

if nargin == 1
    n = 128;
elseif nargin == 2
    n = varargin{1};
else
    error('Too many nargins!')
end

% This next line is a shitshow. I want logarithmically spaced points
% between 1 and the half length of my centres vector. Logspace wants the
% decades (i.e.: 5e5 is in decade 5), so use log10 to get that. Allanvar
% wants a list of unique integers, so round and remove duplicates.
lT = unique(round(logspace(0,floor(log10(nP/2)),n)));
end