function [wb, xb] = blockAverage(w, x, nB, spacing)
% Blocking from Berg-Sørensen 2004 section IV.
if isrow(x) && isvector(x)
    x = x';
end
if isrow(w) && isvector(w)
    w = w';
end

if any(w < w(1))
    error('Please give frequencies in ascending order')
end

if ~exist('spacing','var')
    spacing = 'linear';
end

switch spacing
    case 'linear'
        inds = 1:nB:size(x,1);
    case 'log'
        if nB < 30
            error('Highly recommend nB at least 30, you chose %i', nB)
        end
        inds = [1, 2, 6, 11, 31, 61, 100:100:900, unique(round(logspace(3, log10(size(x,1)), nB-14)))];
        indmax = find(diff(inds)>1e4,1,'first');
        inds = [inds(1:indmax), inds(indmax)+1e4:1e4:size(x,1)];
    case 'log2'
        inds = unique(round([1, 2:10:22, ceil(1.667.^(7:size(x,1)))]));
        inds = inds(inds<size(x,1));
    otherwise
        error('Unrecognised spacing: %s', spacing)
end

if inds(end) ~= size(x,1)
    inds(end) = size(x,1);
end
xb = zeros(length(inds)-1, size(x,2));
wb = zeros(length(inds)-1,1);

% There's a quicker way of doing this using reshape (shrug)
for idx = 1:length(inds)-1
    xb(idx,:) = mean(x(inds(idx):inds(idx+1)-1,:),1);
    wb(idx) = mean(w(inds(idx):inds(idx+1)-1,:),1);
end
