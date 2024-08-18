function [x,t]=diffusive(N,Dt,x1,R,T,eta)
%% [x,t]=diffusive(N,Dt,x1,R,T,eta)
% Brownian dynamics of free diffusing particle in diffusive regime, from
% Volpe 2013 S.I.
%
% Modified by Will Hardiman, August 2024
kB = 1.38e-23;
 % Boltzmann constant [J/K]
gamma = 6*pi*R*eta; % friction coefficient
D = kB*T/gamma;
 % diffusion coefficient
x(1)=x1;
 % initial condition
for i = 2:1:N
x(i) = x(i-1) + sqrt(2*D*Dt)*randn();
end
t = [0:Dt:(N-1)*Dt];

end