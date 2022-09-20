function [z, rho, ngp, N] = accu_Ddist(accumulated, normT)
%% [z, rho, ngp] = accu_Ddist(accumulated, [normT])
% Calculate Ddist for all data in accumulated and return the normalised
% increments, z, the p.d.f., rho, and the non-gaussian parameter α with lag
% times, ngp.

if ~exist('normT', 'var')
    normT = false;
end

% This is a guess at the size we'll need
rho = zeros(140, 28, 2);
% Record the number of samples at each normalised delay
N = zeros(1,size(rho,2)); 
% All scenarios
dt = [1; 2; 4; 6] .* logspace(-4, 3, 8);
dt = dt(:);

% Loop over everything we have here
for dIdx = 1:size(accumulated,2)
    for cIdx = 1:size(accumulated{1,dIdx},2)
        for rIdx = 1:size(accumulated{1,dIdx}{1,cIdx},2)
            % Get the track and calculate the Ddist in a fresh msdanalyzer
            t = accumulated{1,dIdx}{1,cIdx}(rIdx).msd.tracks;
            m = msdanalyzer(1, 'um', 's', 'log');
            m = m.addAll(t);
            m = m.computeDdist;
            
            % Janky but it should work
            if normT
                tnorm = accumulated{1,dIdx}{1,cIdx}(rIdx).tnorm;
            else
                tnorm = 1;
            end
            
            dtnorm = round(m.Ddist{1,1}(:,1),1, 'significant')...
                ./tnorm;
            [tErr, tInd] = min(abs(dt - dtnorm(1)));
            fprintf('tErr = %gs, or %g%% \n', tErr, tErr./dtnorm(1))
            
            % Look I'm really sorry about this, I just couldn't think of a
            % better way to make it robust to unexpected size of Ddist.
            ddsz = size(m.Ddist{1,2}(1:floor(end/2),:))+[0, tInd];
            larger = ddsz > size(rho(:,:,1));
            if any(larger)
                warning('Resizing rho to [%i %i]', ddsz(1), ddsz(2))
            end
            if larger(1)
                rho(ddsz(1),1,1) = 0;
            end
            if larger(2)
                rho(1,ddsz(2),1) = 0;
                N(1,ddsz(2)) = 0;
            end
        
            % Add the counts to the total
            for dim = 1:2
                ind = size(m.Ddist{dim,2},2);
                rho(:,tInd+(1:ind),dim) = rho(:,tInd+(1:ind),dim) + m.Ddist{dim,2}(1:floor(end/2),1:ind);
            end
            % Count how many samples we've added
            try
                N(tInd+(1:ind)) = N(tInd+(1:ind)) + m.Ddist{1}(:,2)';
            catch
                warning('huh')
            end
        end
    end
end

% Normalise and prepare for output
rho = rho ./ sum(rho);
z = m.Ddist{1,2}(ceil(end/2):end,1);
z = z(1:end-1) + diff(z(1:2))/2;

% NGP - see Bursac 2005 Nat.Mat. or Weeks 2002 PRL.
%  <z^4> / (3 <z^2> ^2 ) - 1
alpha = @(z, rho) sum(rho.*(z.^4),1) ./ (3 .* sum(rho.*(z.^2),1) .^2) - 1;
ngp = alpha(z, rho);

try
    ngp = [ngp; repmat(dt(1:size(ngp,2))',1,1,2)];
catch
    warning('Size mismatch in ngp?')
end


