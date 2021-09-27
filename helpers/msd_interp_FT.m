function [omega, G1, G2] = msd_interp_FT(tau, msd, g0, eta, nF, interpF, interpM)
%% Calculate rheological Fourier Transform of msd sampled at times tau.
% Use zero-time point g0, long time inverse gradient eta, calculate at nF
% different frequencies, interpolate msd by factor interpF using method
% interpM.

if ~exist('g0', 'var')
    g0 = 0;
end
if ~exist('ginf', 'var')
    eta = inf;
end
if ~exist('nF', 'var')
    nF = length(tau);
end
if ~exist('interpF','var')
    interpF = 1e3;
end
if ~exist('interpM','var')
    interpM = 'pchip';
end

taui = logspace(log10(tau(1)), log10(tau(end)), length(tau).*interpF)';
msdi = interp1(tau, msd, taui, interpM);

[omega, Gv] = rheoFDFT_Evans_vec(taui, msdi, nF, g0, eta);

G1 = real(Gv);
G2 = imag(Gv);
