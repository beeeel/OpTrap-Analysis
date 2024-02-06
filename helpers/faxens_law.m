function correction = faxens_law(a, h, parallel)
%% correction = faxens_law(a, h, [parallel])
% Calculate viscosity correction based on Faxen's law, default to parallel
% correction = 1 / (1 - 9 / 16 * a / h + 1 / 8 * (a / h)^3 - 45 / 256 * (a / h).^4 - 1 / 16 * (a / h).^5);

if ~exist('parallel','var') || parallel
    correction = 1 ./ (1 - 9 / 16 * a ./ h + 1 / 8 * (a ./ h).^3 - 45 / 256 * (a ./ h).^4 - 1 / 16 * (a ./ h).^5);
else
    correction = 1 ./ (1 - 9 / 8 * a ./ h + 1 / 2 * (a ./ h).^3 );
end
