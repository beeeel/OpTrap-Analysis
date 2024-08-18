function [x,t]=inertial(N,Dt,x1,x2,R,T,eta,d)
%% [x,t]=inertial(N,Dt,x1,x2,R,T,eta,d)
% Brownian dynamics of free diffusing particle in inertial regime, from
% Volpe 2013 S.I.
%
% Modified by Will Hardiman, August 2024
kB = 1.38e-23;
 % Boltzmann constant [J/K]
gamma = 6*pi*R*eta; % friction coefficient
m = 4/3*pi*R^3*d;
 % particle mass
tau = m/gamma;
 % momentum relaxation time
x(1)=x1; x(2)=x2;
 % initial conditions
 
% preallocate
for i = 3:1:N
x(i) = (2+Dt*gamma/m)/(1+Dt*gamma/m)*x(i-1) ...
- 1/(1+Dt*gamma/m)*x(i-2) ...
+ sqrt(2*kB*T*gamma)/(m+Dt*gamma)*randn();
end
t = [0:Dt:(N-1)*Dt];

end