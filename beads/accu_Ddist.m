function [z, rho, ngp, N] = accu_Ddist(accumulated, varargin)
%% [z, rho, ngp] = accu_Ddist(accumulated, [...])
% Calculate Ddist for all data in accumulated and return the normalised
% increments, z, the p.d.f., rho, and the non-gaussian parameter α with lag
% times, ngp.
% 
% Additional options as name-value pairs. Possibilities:
%  normT        Normalise time domain before calculating increments
%  normR        How to normalise spatial domain - either std (Weeks 2002) or bG (me!)
%  posiT        Which times in accumulated to use (-1 for -ve, 0 for all, +1 for +ve)

p = inputParser;

addRequired(p, 'accumulated', @(x) isa(x, 'cell'));
addOptional(p, 'normT', false, @(x) islogical(x) || any(strcmp(x, {'lowT', 'alphaMin', 'highT'})));
addOptional(p, 'normR', 'std', @(x) any(strcmp(x, {'std', 'bG'})));
addOptional(p, 'posiT', 0, @(x) x == -1 || x == 0 || x == 1);
addOptional(p, 'obsTnorm', false, @(x) islogical(x));


parse(p, accumulated, varargin{:});

normT = p.Results.normT;
normR = p.Results.normR;
posiT = p.Results.posiT;
obsTnorm = p.Results.obsTnorm;

% All scenarios
dt = [1; 2; 4; 6] .* logspace(-4, 6, 11);
dt = dt(:)';
% This is a guess at the size we'll need
rho = zeros(140, length(dt), 2);
N = zeros(length(dt),2);
% Count number of nans
nnans = [0 0];
nC = 0;

% Loop over everything we have here
for dIdx = 1:size(accumulated,2)
    for cIdx = 1:size(accumulated{1,dIdx},2)
        for rIdx = 1:size(accumulated{1,dIdx}{1,cIdx},2)
            if obsTnorm
                obstTnf = -mean(accumulated{1,dIdx}{2,cIdx});
            else
                obstTnf = 0;
            end
            obsT = accumulated{1,dIdx}{2,cIdx}(rIdx)+obstTnf;
            if ~posiT || sign(obsT) == posiT
                nC = nC + 1;
                % Get the track and calculate the Ddist in a fresh msdanalyzer
                t = accumulated{1,dIdx}{1,cIdx}(rIdx).msd.tracks;
                m = msdanalyzer(1, 'um', 's', 'log');
                m = m.addAll(t);
                if strcmp(normR, 'std')
                    m = m.computeDdist;
                else
                    normR = kBT(293) ./ accumulated{1,dIdx}{1,cIdx}(rIdx).stiff2(:,1,1);
                    m = m.computeDdist([], normR);
                end
                
                % Janky but it might now work
                if normT
                    if islogical(normT) || strcmp(normT, 'lowT')
                        tnorm = accumulated{1,dIdx}{1,cIdx}(rIdx).tnorm;
                    elseif strcmp(normT, 'highT')
                        tnorm = accumulated{1,dIdx}{1,cIdx}(rIdx).tCH;
                    elseif strcmp(normT, 'alphaMin')
                        msd = cat(3, accumulated{1,dIdx}{1,cIdx}(rIdx).msd.msd{:});
                        msd = msd(1:round(0.8*end), :, :);
                        tau = msd(:,1,1);
                        msd = squeeze(msd(:,2,:));
                        [dydx, tout] = msd_gradientor(tau, msd);
                        [ind, ~] = find(dydx == min(dydx));
                        tnorm = tout(ind)';
                    end
                    
                    if any(isnan(tnorm))
                        warning('found %i NaNs in tnorm for day %i cell %i rep %i', sum(isnan(tnorm)), dIdx, cIdx, rIdx)
                        nnans = nnans + isnan(tnorm);
                    end
                else
                    tnorm = [1 1];
                end
                
                dtnorm = round(m.Ddist{1,1}(:,1),1, 'significant')...
                    ./tnorm;
                [tErr, tInd] = min(abs(dt - dtnorm(1,:)'), [], 2);
                [~, ind] = max(tErr);
                fprintf('tErr = %gs, or %g%% \n', tErr(ind), tErr(ind)./dtnorm(1,ind))
                
                % Look I'm really sorry about this, I just couldn't think of a
                % better way to make it robust to unexpected size of Ddist.
                ddsz = size(m.Ddist{1,2}(1:floor(end/2),:))+[0, max(tInd)];
                larger = ddsz > size(rho(:,:,1));
                if any(larger)
                    warning('Resizing rho to [%i %i]', ddsz(1), ddsz(2))
                end
                if larger(1)
                    rho(ddsz(1),1,1) = 0;
                end
                if larger(2)
                    rho(1,ddsz(2),1) = 0;
                    N(ddsz(2), 1) = 0;
                end
                
                try 
                    % Add the counts to the total
                    for dim = 1:2
                        if ~isnan(tnorm(dim))
                            ind = size(m.Ddist{dim,2},2);
                            N(tInd(dim)+(1:ind),dim) = N(tInd(dim)+(1:ind),dim) + m.Ddist{1}(:,2);
                            rho(:,tInd(dim)+(1:ind),dim) = rho(:,tInd(dim)+(1:ind),dim) + m.Ddist{dim,2}(1:floor(end/2),1:ind);
                        end
                    end
%                     disp('dt and N this time:')
%                     str = repmat('\n%.1g\t\t%g\t\t%g\n',1,size(N,1));
%                     a = zeros(2,size(N,1));
%                     a(:,tInd(dim)+(1:ind)) = round(m.Ddist{1}(:,1:2)',1,'significant');
%                     a = [a; N(:,1)'];
%                     fprintf(str, a) 
%                     fprintf('dt\t\tN now \t\t N total\n\n')
%                     if N(19,1)==0
%                         error('ya')
%                     end
                catch ME 
                    error(ME.message, 'uhhh')
                end
            end
        end
    end
end

% Normalise and prepare for output
rho = rho ./ reshape(N,1,[],2);
z = m.Ddist{1,2}(ceil(end/2):end,1);
z = z(1:end-1) + diff(z(1:2))/2;

if any(nnans > 0)
    fprintf('\n\n\t\tFailed time normalisation for [%i, %i] out of %i in respective dimensions\n\n', nnans, nC)
end

% NGP - see Bursac 2005 Nat.Mat. or Weeks 2002 PRL.
%  <z^4> / (3 <z^2> ^2 ) - 1
alpha = @(z, rho) sum(rho.*(z.^4),1) ./ (3 .* sum(rho.*(z.^2),1) .^2) - 1;
ngp = alpha(z, rho);

try
    ngp = [ngp; repmat(dt,1,1,2)];
catch
    warning('Size mismatch in ngp?')
end


