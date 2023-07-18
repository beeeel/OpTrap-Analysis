function [M, varargout] = mmsd(varargin)
%% [M, [kappa]] = mmsd(msdObj, tR)
%% [M, [kappa]] = mmsd(tau, msd, tR)
% Calculate mean MSD and optionally stiffness (κ = 2kBT/MSD) given a time
% range to average over. Either give an msdanalyzer or tau and msd. In the
% latter case, multiple observations should go in third dimension; tau can
% either be reused (in which case size(tau,3) == 1) or one tau for each msd
% (in which case size(tau,3) == size(msd,3)).
%
% If you want kappa to have sensible answers, give msd in m^2, or multiply
% kappa by 10^12.

if nargin == 2
    msdObj = varargin{1};
    if ~isa(msdObj, 'msdanalyzer')
        error('%i input arguments given, expected first to be msdanalyzer, instead got %s', nargin, class(msdObj))
    end
    tR = varargin{2};
    tmp = cat(3, msdObj.msd{:});
    tau = tmp(:,1,:);
    msd = tmp(:,2,:);
    
    if strcmp(msdObj.space_units, 'um')
        mult = 1e-12;
    end
else
    tau = varargin{1};
    msd = varargin{2};
    tR = varargin{3};
    mult = 1;
    if size(tau,2) ~= 1 || size(msd,2) ~= 1 || (size(msd,3) ~= size(tau,3) && size(tau,3) ~= 1)
        error('Check input sizes for tau and msd. Should be 1 row per delay time, 1 column, and put repeats in the third dimension')
    end
end

M = zeros(1,size(msd,3));

if numel(tR) == 2
    inds = tau > tR(1) & tau < tR(2);
elseif size(tR,3) == size(tau,3)
    inds = false(size(msd));
    for idx = 1:size(tau,3)
        inds(:,1,idx) = tau(:,1,idx) > tR(1) & tau(:,1,idx) < tR(2);
    end
elseif size(tR,3) == size(msd,3)
    inds = false(size(msd));
    for idx = 1:size(msd,3)
        inds(:,1,idx) = tau > tR(1) & tau < tR(2);
    end
else 
    error('Check your tR sizes compared to msd and tau sizes')
end

for idx = 1:size(msd,3)
    if size(tau,3) == size(msd,3)
        M(idx) = mean(msd(inds(:,:,idx),:,idx),1);
    else
        M(idx) = mean(msd(inds),1);
    end
end

if nargout == 2
    varargout{1} = 2 * kBT(298) ./ (mult * M);
end