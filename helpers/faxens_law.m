function correction = faxens_law(a, h)
%% correction = faxens_law(a, h)
% Calculate viscosity correction based on Faxen's law
% correction = 1 / (1 - 9 / 16 * a / h + 1 / 8 * (a / h)..^3 - 45 / 256 * (a / h).^3 - 1 / 16 * (a / h).^5);

correction = 1 ./ (1 - 9 / 16 * a ./ h + 1 / 8 * (a ./ h).^3 - 45 / 256 * (a ./ h).^3 - 1 / 16 * (a ./ h).^5);
