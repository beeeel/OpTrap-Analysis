function [omega, G1, G2, varargout] = msd_interp_FT(tau, msd, g0, eta, nF, interpF, interpM, interpSpacing)
%% Calculate rheological Fourier Transform of msd sampled at times tau.
% [omega, G1, G2, [MFT1, MFT2]] = msd_interp_FT(tau, msd, g0, eta, nF, interpF, interpM, interpSpacing)
% Use zero-time point g0, long time inverse gradient eta, calculate at nF
% different frequencies, interpolate msd by factor interpF using method
% interpM. Space interpolated points log or linear.

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
    interpM = 'spline';
end
if ~exist('interpSpacing','var')
    interpSpacing = 'log';
    allow0FT = false;
else
    allow0FT = true;
    if any(tau == 0)
        msd = msd(tau ~= 0);
        tau = tau(tau ~= 0);
    end
end



switch interpSpacing
    case 'log'
        taui = logspace(log10(tau(1)), log10(tau(end)), length(tau).*interpF)';
    case 'linear'
        taui = linspace(tau(1), tau(end), length(tau).*interpF)';
%         allow0FT = true;
    otherwise
        error('Unrecognised interpSpacing: %s', interpSpacing)
end
msdi = interp1(tau, msd, taui, interpM);

[omega, Gv, MFT] = rheoFDFT_Evans_vec(taui, msdi, nF, g0, eta, allow0FT);

G1 = real(Gv);
G2 = imag(Gv);

if nargout >= 5
    varargout = {real(MFT), imag(MFT)};
end