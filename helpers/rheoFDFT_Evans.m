function [varargout] = rheoFDFT_Evans(tau, msd, nOmegas, J0, eta)
%[omega, Y, [Z]] = rheoFDFT_Evans(tau, msd, nOmegas, J0, eta)
%% Compute finite discrete fourier transform using Evans et al 2009
% Requires consideration of the long-time limit of the gradient = 1/η
% Needs to be supplied the zero-time limit of the compliance/MSD, J0

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

% Setup: Frequency sampling
omega = 2*pi./logspace(log10(tau(1)), log10(tau(end)), nOmegas)';
% omega = 2*pi./linspace(tau(1), tau(end), nOmegas);
% omega = 2*pi./tau;

%% Calculate the Fourier transform using eq. 10
Y = 1i * omega * J0 ... % First term
    + ( 1 - exp(-1i * omega * tau(1)) ) * ( msd(1) - J0 ) / tau(1) ... % Second term
    + exp(-1i * omega * tau(end)) / eta; % Third term

% Summation term - calculate summation in vectorized form, for each
% frequency sequentially (embarassingly parallelizable)
for wIdx = 1:nOmegas
    Y(wIdx) = Y(wIdx) + sum( ( diff( msd ) ./ diff( tau ) ) ...
        .* ( exp( -1i * omega(wIdx) * tau( 1:N-1 )) - exp( -1i * omega(wIdx) * tau( 2:N )) ));
end

Z = Y./(-omega.^2); % The Fourier transform of MSD
Y = 1i .* omega ./ Y;

varargout = {omega, Y};
if nargout == 3
    varargout{3} = Z;
end
fprintf('\rDone!%s\n',repmat(' ', 50,1))

end