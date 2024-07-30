function [varargout] = rheoFDFT_Evans_vec(tau, msd, nOmegas, J0, eta, allowZero)
% Y = rheoFDFT_Evans_vec(tau, msd, nOmegas, J0, eta, allowZero)
% [omega, [Y, [Z]]] = rheoFDFT_Evans_vec(tau, msd, nOmegas, J0, eta, allowZero)
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
if ( tau(1) == 0 || min(msd(1,:) == 0) ) && (~exist('allowZero','var') || ~allowZero)
    error('Expected non-zero lag time and MSDs')
end

% Setup: Frequency sampling
if (~exist('allowZero','var') || ~allowZero)
    omega = 1./logspace(log10(tau(1)), log10(tau(end)), nOmegas)';
else
    omega = 2*pi./tau;
    if any(tau == 0)
        omega = omega(tau ~= 0);
        tau = tau(tau ~= 0);
        N = length(tau);
    end
    omega = omega(unique(round(linspace(1, length(omega), nOmegas))));
end

%% Find out how much this can be parallelized
% Determine max array size by asking for an array of increasing size
try
    % First see if we can do it fully vectorized (most memory intense)
    [~] = zeros(2*N, nOmegas);
    % Indexes for summation term
    idx = 1;
    edx = N;
    nEl = N;
catch
    % If that fails, incrementally ask for bigger arrays until it fails
    n = 0;
    m = 1;
    state = true;
    while state && (m*10^n) < nOmegas^2
        % Big steps
        n = n + 1;
        try
            [~] = zeros(nOmegas, m*10^n);
        catch
            n = n - 1;
            while state
                % Small steps
                m = m + 0.1;
                try
                    [~] = zeros(nOmegas, m*10^n);
                catch
                    m = m - 0.1;
                    state = false;
                end
            end
            
        end
    end
    % Final number of elements we can ask for
    nEl = floor(m * 10^n);
    fprintf('Parallelizing factor %i, largest array needed %gGB\n', nEl, nEl*nOmegas*4/1e9);
    
    % Indexes for summation term
    idx = 1;
    if idx + nEl < N
        edx = idx + nEl;
    else
        edx = N;
    end
end
%% Calculate the Fourier transform using eq. 10
Y = 1i * omega * J0 ... % First term
    + ( 1 - exp(-1i * omega * tau(1)) ) * ( msd(1) - J0 ) / tau(1) ... % Second term
    + exp(-1i * omega * tau(end)) / eta; % Third term

% Summation term
while idx < N
    
    fprintf('\r %i:%i / %i', idx, edx, N)
    Y = Y + sum( diff(msd(idx:edx))' ...
        .* ( exp(-1i * omega .* tau(idx:edx - 1)') - exp(-1i * omega .* tau(idx+1:edx)') ) ...
        ./ ( diff(tau(idx:edx))' ), 2);
    idx = edx;
    edx = min(idx + nEl, N);
end


Z = Y./(-omega.^2); % The Fourier transform of MSD
Y = 1 ./ ( 1i * omega .* Z ); % The complex modulus encoded in MSD

if nargout == 2
    varargout = {omega, Y};
elseif nargout == 3
    varargout = {omega, Y, Z};
else
    error('Wrong number of nargouts: %i. Should be 2 or 3', nargout)
end

fprintf('\rDone!\n')

end
