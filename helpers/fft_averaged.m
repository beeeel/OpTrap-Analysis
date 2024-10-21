function [X, w] = fft_averaged(x, dt, avgMode, nAvg, ZP)
%%[X, w] = fft_averaged(x, dt, avg_mode, nAvg, [ZP])
% Calculate fft by subsampling data in time domain, performing fft on each
% shorter series, and averaging the spectra.
%
% 
% Average modes:
%   block: 
%       Bartlett's method

%% Input handling
if ~isvector(x) || ~isscalar(dt)
    error(['x needs to be a vector, size(x) = [%s],' ...
        'and t needs to be scalar, size(t) = [%s]'], ...
        num2str(size(x)), num2str(size(dt)))
end

nT = length(x);
nX = nT/nAvg;

if nX ~= round(nX)
    error('nAvg (%i) needs to be a factor of number of points (%i)', nAvg, nT)
end

%% Reshape data to put new series in the columns of x
if ~iscolumn(x)
    x = x(:);
end

switch avgMode
    case 'Bartlett'
        % Take each block of nT/nAvg points as shorter series.
        % Preserves highest frequencies
        x = reshape(x, nT/nAvg, []);
        X = fft(x,ZP,1);

        X = mean(abs(X(1:end/2,:)).^2,2)/nX;
    case 'subsample'
        % Take one point in every nAvg as shorter series
        % Preserves lowest frequencies
        x = reshape(x, nAvg, [])';
        % Need to resample t so the frequency calculates correctly
        dt = dt*nAvg;
    otherwise
        error('ya goof')
end

if ~exist('ZP','var')
    ZP = size(x,1);
end

%% Do the FFT and average

Fs = 1/dt;
w = linspace(0,Fs/2, ZP/2)';


