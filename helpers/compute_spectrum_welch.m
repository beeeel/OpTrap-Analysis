function [w, X] = compute_spectrum_welch(x, t, segment_length, overlap, use_dpss, alpha)
%% [w, X] = compute_spectrum_welch(x, t, segment_length, overlap, use_dpss, alpha)
% Compute the PSD using Welch's method with the option to use a DPSS window
%
% Inputs:
%   x              - Time-domain signal (vector)
%   t              - Time vector corresponding to the signal x (vector)
%   segment_length - Length of each segment for Welch's method (integer)
%   overlap        - Number of points to overlap between segments (integer)
%   use_dpss       - (Optional) Boolean flag to use DPSS window (true/false)
%   alpha          - (Optional) Time-bandwidth product (used if DPSS window is selected)
%
% Outputs:
%   w - Frequency vector (in Hz)
%   X - Power spectral density (PSD) estimate

% Ensure x and t are row vectors
x = x(:)';
t = t(:)';

% Sampling interval and frequency
dt = t(2) - t(1);   % Assume uniform sampling
fs = 1 / dt;        % Sampling frequency

% Set default values for optional inputs
if nargin < 5 || isempty(use_dpss)
    use_dpss = false;  % Default to not using DPSS window
end

if use_dpss
    if nargin < 6 || isempty(alpha)
        error('When using DPSS window, you must provide the alpha value (time-bandwidth product).');
    elseif alpha < 0 || alpha >= segment_length/2
        error('alpha must be in range [0, %i), you gave %i', segment_length/2, alpha)
    end
    % Use DPSS window with given time-bandwidth product (alpha)
    M = segment_length;  % DPSS window length
    [window, ~] = dpss(M, alpha);
else
    % Use a default window (Hamming)
    window = hamming(segment_length);
end

% Compute the PSD using Welch's method
[X, w] = pwelch(x, window, overlap, [], fs, 'onesided');

% X is the PSD estimate, and w is the corresponding frequency vector
end
