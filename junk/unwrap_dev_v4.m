% DEGREES!!!!! 

% This is the derivative of Gaussian equation with the radial equation for
% an ellipse edge substituded for mu, the centring parameter.
FitEqn = @(a, b, phi, sigma, amp, x, y) ...
                amp .* (y - sqrt((a .* cos(0.0175.*x + phi)).^2 + (b .* sin(0.0175.*x + phi)).^2)) .* ...
                exp(-((y - sqrt((a .* cos(0.0175.*x + phi)).^2 + (b .* sin(0.0175.*x + phi)).^2)).^2)./(2 .* sigma .^2)) ./ ...
                (sqrt(2 .* pi) * sigma.^3);

r = 1:size(Unwrapped,1);
theta = 1:360;
a = 100;
b = 80;
phi = 0;
sigma = 7;
amp = 5e4;

Simulated = FitEqn(a, b, phi, sigma, amp, theta, r');
imagesc(Simulated)
%%
LB = [0, 0, -pi/2, 0, 0];
UB = [200, 200, pi/2, 100, inf];
Start = [100, 100, 0, 5, 1e4];
[ThOut, rOut, IOut] = prepareSurfaceData(theta, r, Unwrapped(:,:,1));
FitObj = fit([ThOut, rOut], IOut,FitEqn,'Upper',UB,'Lower',LB,'Start',Start);
plot(FitObj, [ThOut, rOut], IOut)