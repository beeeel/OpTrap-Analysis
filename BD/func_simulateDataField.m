function [varargout] = func_simulateDataField(opts)
%% [tracks, stiffnesses, msdanalyzer] = func_simulateDataField(opts)
%% Perform BD of particle in optical trap and electric field
% Input an opts structure from func_BDopts, outputs optional - tracks for
% msdanalyzer, stiffnesses (kBT/variance) as [kx1 ky1 kz1; kx2...],
% msdanalyzer with msd calculated.

Nb = opts.Nbeads;

r = opts.radius;                       % radius of bead
t = opts.temp;                         % temperature in Celsius
T = 273+t;                      % temperature in Kelvin

kx = opts.kappaNm(:,1);                      % trap stiffness along x
ky = opts.kappaNm(:,2);                      % trap stiffness along y
kz = opts.kappaNm(:,3);                      % trap stiffness along y


E = opts.E_func;		% Force field defined via function of position
q = 1.6e-19 * opts.q_bead;			% Charge of particle in Coulombs

eta = opts.eta;			% Viscosity of solution

%viscosity-temperature eq for water 
%eta = 2.414E-5 * 10.^(247.8 ./(T -140)); 
%complex viscosity = dynamic viscosity + out-of-phase viscosity (for
%water = dynamic viscosity)

N = opts.Nt;                       % number of samples
fps = 1/opts.dt;                      % sample rate
gamma0 = 6*pi*r.*eta;            % viscous drag force
deltat = 1./fps;                 % duration of each frame

kB=1.38E-23;                    % Boltzmann's constant
D=(kB.*T)./gamma0;                  % diffusion coefficient

% initialise bead position vectors
pos0 = opts.pos0;
x = zeros(Nb,N)+pos0(:,1);
y = zeros(Nb,N)+pos0(:,2);
z = zeros(Nb,N)+pos0(:,3);
fs = zeros(3,N);
time = (0:N-1).*opts.dt;

% random Gaussian noise
rng(opts.rng_seed);
noise = randn(Nb, 3, N);

%% Simulation


for i=1:N-1
    
    
    Ef = E(x(:,i), y(:,i), z(:,i), time(i)) .* q;
    % Calculate position (Langevin's equation) - Volpe and Volpe 2012, 2014
    % restoring force
    x(:,i+1)= x(:,i) -(kx.*deltat.*x(:,i)./gamma0) ...	% Trap force
        + noise(:,1,i).*sqrt(2*D*deltat) ... 			% Diffusion
        + Ef(:,1)*deltat./gamma0;                         % Electric force

    fs(:,i+1) = [-(kx(1).*deltat.*x(1,i)./gamma0(1)), ...	% Trap force
        + noise(1,1,i).*sqrt(2*D(1)*deltat), ... 			% Diffusion
        + Ef(1,1)*deltat./gamma0(1)];
    
    y(:,i+1)= y(:,i) -(ky.*deltat.*y(:,i)./gamma0) ...
        + noise(:,2,i).*sqrt(2*D*deltat) ...
        + Ef(:,2)*deltat./gamma0;
    
    z(:,i+1)= z(:,i) -(kz.*deltat.*z(:,i)./gamma0) ...
        + noise(:,3,i).*sqrt(2*D*deltat) ...
        + Ef(:,3)*deltat./gamma0;
    
end

% kxe(j) = kB*T/var(x(:,j));
% kye(j) = kB*T/var(y(:,j));


%% 
tracks = cell(1,Nb);
for j = 1:Nb
	tracks{j} = [time',x(j,:)',y(j,:)', z(j,:)'];
end
varargout{1} = tracks;

if nargout > 1
    kxe = kB*T./var(x,0,2);
    kye = kB*T./var(y,0,2);
    kze = kB*T./var(z,0,2);
    varargout{2} = [kxe kye kze];
end

if nargout > 2
    msd = msdanalyzer(3, 'm','s','log');
    
    msd = msd.addAll(tracks);
    msd = msd.computeMSD;
    varargout{3} = msd;
end
