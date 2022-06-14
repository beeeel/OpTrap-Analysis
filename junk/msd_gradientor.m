function [dydx, varargout] = msd_gradientor(tau, msd, varargin)
%% Calculates gradients in log-log space
% [dydx, tout] = msd_gradientor(tau, msd, [method, nP])
% Can use 2-point average, linear fitting or least squares from Ling 2019

if isempty(varargin) || isempty(varargin{1})
    method = 'lsq';
else
    method = varargin{1};
end
if nargin == 4
    nP = varargin{2};
end

switch method
    case 'piecewise'
        % Calculate pointwise gradient in log space
        tRs = tau(2:end)./tau(1:end-1);
        mRs = msd(2:end,:)./msd(1:end-1,:);
        
        dydx = log(mRs) ./ log(tRs);
        
        % Rolling average
        if ~exist('nP','var')
            nP = 15;
        end
        kern = ones(nP,1);
        dydx = conv(dydx, kern, 'valid')./nP;
        tout = conv(tau(1:end-1), kern, 'valid') ./ nP;
    case 'fitting'
        if ~exist('nP','var')
            nP = 40;
        end
        dydx = zeros(size(msd,1)-nP, size(msd,2));
        tout = nan(length(tau)-nP, 1);
        tau = tau(msd > 0); % Careful, if MSD is a matrix this will shit itself
        msd = msd(msd > 0);
        for jdx = 1:size(msd,2)
            for idx = 1:length(tau)-nP
                tdata = tau(idx:idx+nP);
                mdata = msd(idx:idx+nP,jdx);
                
                tout(idx) = mean(tdata);
                
                fo = fit(log(tdata), log(mdata), 'Poly1');
                dydx(idx, jdx) = fo.p1;
            end
        end
    case 'lsq'
        if ~exist('nP','var')
            nP = 15;
        end
        msd = reshape(msd(tau > 0,:), [], size(msd,2));
        tau = tau(tau > 0); % Careful, if MSD is a matrix this will shit itself (well, hopefully not any more)
        assert(size(tau,1) == size(msd,1), 'Ugh you silly, tau and MSD need to have same number of rows')
        dydx = zeros(size(msd,1)-nP, size(msd,2));
        tout = nan(length(tau)-nP, 1);
        for jdx = 1:size(msd,2)
            for idx = 1:length(tau)-nP
                tdata = tau(idx:idx+nP);
                mdata = msd(idx:idx+nP,jdx);
                
                tout(idx) = mean(tdata);
                
                [alpha, ~] = leastSq(log(tdata), log(mdata));
                dydx(idx, jdx) = alpha;
            end
        end
        
    otherwise
        error('Choose a correct method')
end

if nargout == 2
    varargout = {tout};
elseif nargout > 2
    error('Too many nargouts in msd_gradientor')
end

end