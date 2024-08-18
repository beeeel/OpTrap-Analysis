function [x,y,z,t]=trapped(N,Dt,x1,y1,z1,R,T,eta,kx,ky,kz)
%% [x,y,z,t]=trapped(N,Dt,x1,y1,z1,R,T,eta,kx,ky,kz)
% Brownian dynamics of optically trapped particle in diffusive regime, from
% Volpe 2013 S.I.
%
% Modified by Will Hardiman, August 2024
kB = 1.38e-23;
 % Boltzmann constant [J/K]
gamma = 6*pi*R*eta; % friction coefficient
D = kB*T/gamma;
 % diffusion coefficient
x(1)=x1;y(1)=y1;z(1)=z1;
 % initial condition
for i = 2:1:N
% Deterministic step
x(i) = x(i-1) - kx*Dt/gamma*x(i-1);
y(i) = y(i-1) - ky*Dt/gamma*y(i-1);
z(i) = z(i-1) - kz*Dt/gamma*z(i-1);
% Diffusive step
x(i) = x(i) + sqrt(2*D*Dt)*randn();
y(i) = y(i) + sqrt(2*D*Dt)*randn();
z(i) = z(i) + sqrt(2*D*Dt)*randn();
end
t = [0:Dt:(N-1)*Dt];

end