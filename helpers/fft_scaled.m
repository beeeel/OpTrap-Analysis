function [w, X] = fft_scaled(t, x, varargin)
%% [w, X] = fft_scaled(t, x, varargin)
% Calculate Fourier transform of xf, return one-sided spectra in m and
% frequency vector in Hz. Assumes x in m and t in s. Plots X in m and w in
% Hz.

doPlots = true;
fn = '';

if nargin >= 3
    doPlots = varargin{1};
end
if nargin >= 4
    ax = varargin{2};
else
    figure
    ax = axes;
end
if nargin >= 5
    fn = varargin{3};
end

n_points = length(x);

% Get the FFT
X = fft(x);
X = abs(eval([fn '(X/n_points)']));
X = X(1:floor(end/2)+1);
X(2:end-1) = 2*X(2:end-1);

% Frequency in index (k+1) is k cycles per whole dataset, so convert to Hz
w = (0:ceil((n_points-1)/2))./diff(t([1 end]));

if doPlots
    plot(ax, w, X)
    
    title('Frequency spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (m)')
    
    drawnow
end