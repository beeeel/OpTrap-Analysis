function [w, X] = fft_scaled(t, x, varargin)
%% [w, X] = fft_scaled(t, x, [doPlots, ax, fn, zp])
% Calculate Fourier transform of xf, return one-sided spectra in m and
% frequency vector in Hz. Assumes x in m and t in s. [Optional] Plots X in
% m and w in Hz, [optional] on axis ax, and evaluate function fn using
% input X. [Optional] apply zero padding when doing the FFT

doPlots = true;
fn = '';

if nargin >= 3
    doPlots = varargin{1};
end
if nargin >= 4 && ~isempty(varargin{2})
    ax = varargin{2};
elseif doPlots
    figure
    ax = axes;
end
if nargin >= 5 && ~isempty(varargin{3})
    fn = varargin{3};
end
if nargin >= 6 && ~isempty(varargin{4})
    zp = varargin{4};
    if ~isfinite(zp) || zp < 0
        warning('Overriding zero padding in fft_scaled')
        zp = [];
    end
else 
    zp = [];
end

n_points = size(x,2);

if n_points == 1
    n_points = size(x,1);
    x = x.';
end


if zp < n_points
    warning('Overriding zero padding option: %g', zp)
    zp = n_points;
elseif isempty(zp)
    zp = n_points;
else
    warning('Zero padding was seen causing an error of ~ a factor of 3')
    tmp = zeros(1,2^nextpow2(zp));
    inds = [zp/2-ceil(n_points/2), zp/2+floor(n_points/2)];
    tmp(inds(1):inds(2)-1) = x;
    tmp(1:inds(1)-1) = x(1);
    tmp(inds(2):end) = x(end);
    x = tmp;
    clear tmp
end

%% Get the FFT
X = fft(x, zp, 2);
X = eval([fn '(X/zp)']);
X = X(:,1:floor(end/2)+1);
X(:,2:end-1) = 2*X(:,2:end-1);

% Sampling frequency
Fs = 1./diff(t([1 2]));
w = 0:Fs/zp:Fs/2;

if doPlots
    plot(ax, w, abs(X))
    
    title('Frequency spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (m)')
    
    drawnow
    set(ax,'XScale','log','YScale','log')
end
