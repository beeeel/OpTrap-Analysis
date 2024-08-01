function [wb, xb] = blockAverage(w, x, nB)
% Blocking from Berg-Sørensen 2004 section IV.
if isrow(x)
    x = x';
end
if isrow(w)
    w = w';
end

inds = 1:nB:size(x,1);
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
