function [Y] = rheoFDFT_Evans_vec(tau, msd, J0, eta)
%Y = rheoFDFT_Evans_vec(tau, msd, J0, eta)
%% Compute finite discrete fourier transform using Evans et al 2009
% Requires consideration of the long-time limit of the gradient = 1/η
% Needs to be supplied the zero-time limit of the compliance/MSD, J0
%
% Do calculation in vectorized manner

% Input checking
N = length(tau);
try
    validateattributes(msd, {'numeric'}, {'nrows', N})
catch
    error('Received msd with %i rows, and tau with %i rows', size(msd,1), size(tau,1))
end

% If first tau or msd is zero, throw an error
if tau(1) == 0 || min(msd(1,:) == 0)
    error('Expected non-zero lag time and MSDs')
end

% Setup
omega = 2*pi./tau;

%% Find out how much this can be parallelized
% Determine max array size by asking for an array of increasing size
n = 0;
m = 1;
state = true;
while state && (m*10^n) < N^2
    n = n + 1;
    try
        [~] = zeros(N, m*10^n);
    catch
        n = n - 1;
        while state
            m = m + 0.1;
            try
                [~] = zeros(N, m*10^n);
            catch
                m = m - 0.1;
                state = false;
            end
        end
        
    end
end
nEl = floor(m * 10^n);
fprintf('Parallelizing factor %i, largest array needed %gGB\n', nEl, nEl*N*4/1e9);
%% Calculate the Fourier transform using eq. 10
Y = 1i * omega * J0 ... % First term
    + ( 1 - exp(-1i * omega * tau(1)) ) * ( msd(1) - J0 ) / tau(1) ... % Second term
    + exp(-1i * omega * tau(end)) / eta; % Third term
% Summation term
idx = 1;
if idx + nEl < N
        edx = idx + nEl;
    else
        edx = N;
end
while idx <= N
    
    fprintf('\r %i / %i', idx, N)
    Y = Y + sum(( diff(msd(idx:edx))' ) ...
        .* ( exp(-1i * omega .* tau(idx:edx - 1)') - exp(-1i * omega .* tau(idx+1:edx)') ) ...
        ./ ( diff(tau(idx:edx))' ), 2);
    idx = edx + 1;
    edx = min(idx + nEl, N);
end


Y = 1i .* omega ./ Y;

fprintf('\rDone!\n')

end
