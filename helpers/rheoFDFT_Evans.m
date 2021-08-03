function [Y] = rheoFDFT_Evans(tau, msd, J0, eta)
%Y = rheoFDFT_Evans(tau, msd, omega)
% Compute finite discrete fourier transform using Evans et al 2009
%% Requires consideration of the long-time limit of the gradient = 1/η

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

%% Calculate the Fourier transform using eq. 10
Y = 1i * omega * J0 ... % First term
    + ( 1 - exp(-1i * omega * tau(1)) ) * ( msd(1) - J0 ) / tau(1) ... % Second term
    + exp(-1i * omega * tau(end)) / eta; % Third term
% Summation term
for k = 2:N
    fprintf('\r%i / %i', k, N)
    Y = Y + ( msd(k) - msd(k - 1) ) ...
        .* ( exp(-1i * omega * tau(k - 1)) - exp(-1i * omega * tau(k)) ) ...
        ./ ( tau(k) - tau(k - 1) );
end

Y = 1i .* omega ./ Y;

fprintf('\rDone!\n')

end