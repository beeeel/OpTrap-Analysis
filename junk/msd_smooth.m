function Jt = msd_smooth(t, eta, varargin)
%% Jt = msd_smooth(t, eta, [J1, tau1, J2, tau2, ... ])
% Compliance/MSD function as a sum of exponentials. Number of exponentials
% determined by number of Js and taus (need equal numbers of each)
% See "Analysis of the linear viscoelasticity of polyelectrolytes by
% magnetic microrheometry", Tassieri et al (2010), eq. 11.

if mod(length(varargin), 2)
    error('Need an even number of Js and taus')
end

Jt = t./eta;
for n = 1:( length(varargin) / 2 )
    Jt = Jt + varargin{2*n-1} .* (1 - exp(-t./varargin{2*n}));
end

end