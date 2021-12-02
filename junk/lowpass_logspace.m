function G = lowpass_logspace(omega, G, lpFrq, varargin)
%% G = lowpass_logspace(omega, G, lowpassFrq, [nInterp])

% nInterp = 1e4;
% if nargin > 3
%     nInterp = varargin{1};
% end
    
% omega1 = linspace(omega(1), omega(end), nInterp);
% G1 = interp1(omega, G, omega1);

% G1 = lowpass(G1, lpFrq, 1/diff(omega1([2 1])));
G = lowpass(G, lpFrq);

% G = interp1(omega1, G1, omega);
