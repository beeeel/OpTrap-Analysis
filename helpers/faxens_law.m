function correction = faxens_law(a, h, parallel)
%% correction = faxens_law(a, h, [parallel])
% Calculate viscosity correction based on Faxen's law, default to parallel
% Uses corrections from Schaffer 2007:
% parallel correction = 1 / (1 - 9 / 16 * a / h + 1 / 8 * (a / h)^3 - 45 / 256 * (a / h).^4 - 1 / 16 * (a / h).^5);
% perpendicular = 1 / (1 - 9/8 * a/h + 1/2 * (a/h)^3 - 57/100 * (a/h)^4 + 1/5 * (a/h)^5 + 7/200 * (a/h)^11 - 1/25 * (a/h)^12)

if ~exist('parallel','var') || parallel
    correction = 1 ./ (1 - 9 / 16 * a ./ h + 1 / 8 * (a ./ h).^3 - 45 / 256 * (a ./ h).^4 - 1 / 16 * (a ./ h).^5);
else
    % V1: Not sure where this is from
    % correction = 1 ./ (1 - 9 / 8 * a ./ h + 1 / 2 * (a ./ h).^3 );

    % V2: From Schaffer 2007
    correction = 1 ./ (1 - 9/8 * a./h + 1/2 * (a./h).^3 - 57/100 * (a./h).^4 + 1/5 * (a./h).^5 + 7/200 * (a./h).^11 - 1/25 * (a./h).^12);
end
