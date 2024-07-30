function [varargout] = func_simulateDataField_v2(opts)
%% [tracks, Xforces, stiffnesses, msdanalyzer] = func_simulateDataField_v2(opts)
%% Perform BD of particle in optical trap and surface force, use Faxen's law for viscosity
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
% q = 1.6e-19 * opts.q_bead;			% Charge of particle in Coulombs

eta = opts.eta;			% Viscosity of solution

%viscosity-temperature eq for water 
%eta = 2.414E-5 * 10.^(247.8 ./(T -140)); 
%complex viscosity = dynamic viscosity + out-of-phase viscosity (for
%water = dynamic viscosity)

N = opts.Nt;                       % number of samples
deltat = opts.dt;                 % duration of each frame
fps = 1/deltat;                      % sample rate


kB=1.38E-23;                    % Boltzmann's constant
gamma0 = 6*pi*r.*eta;            % viscous drag force
fax = @(h) [1 1 0] .* faxens_law(r, h, true) + [0 0 1] .* faxens_law(r, h, false);
% Check timestep vs trap stiffness
if any(gamma0./opts.kappaNm < 10/fps)
    warning('Large time step compared to trap relaxation time')
end
% initialise bead position vectors
pos0 = opts.pos0;
x = zeros(Nb,N)+pos0(:,1);
y = zeros(Nb,N)+pos0(:,2);
z = zeros(Nb,N)+max(pos0(:,3), r*1.001);
fs = zeros(3, N);
Efs = zeros(3, N);
time = (0:N-1).*opts.dt;

% random Gaussian noise
if opts.rng_seed ~= 0; rng(opts.rng_seed); end
noise = randn(Nb, 3, N);

%% Simulation


for i=1:N-1

    gamma0 = 6*pi*r.*eta.*fax(z(:,i));            % viscous drag force
    D=(kB.*T)./gamma0;                  % diffusion coefficient

    Ef = E(x(:,i), y(:,i), z(:,i), time(i));
    % Calculate position (Langevin's equation) - Volpe and Volpe 2012, 2014
    % restoring force
    x(:,i+1)= x(:,i) - (kx.*(x(:,i)-pos0(:,1)).*deltat./gamma0(1)) ...	% Trap force
        + noise(:,1,i).*sqrt(2*D(1)*deltat) ... 			% Diffusion
        + Ef(:,1)*deltat./gamma0(1);                         % Electric force
    
    y(:,i+1)= y(:,i) -(ky.*deltat.*(y(:,i)-pos0(:,2))./gamma0(2)) ...
        + noise(:,2,i).*sqrt(2*D(2)*deltat) ...
        + Ef(:,2)*deltat./gamma0(2);
    
    z(:,i+1)= z(:,i) -(kz.*deltat.*(z(:,i)-pos0(:,3))./gamma0(3)) ...
        + noise(:,3,i).*sqrt(2*D(3)*deltat) ...
        + Ef(:,3)*deltat./gamma0(3);
    
    fs(:,i) = [-(kz(1).*(z(1,i)-pos0(:,3))), ...	% Trap force
        noise(1,3,i).*sqrt(2*D(3)*deltat).*gamma0(3)/deltat, ... 			% Diffusion
        Ef(1,3)];
    Efs(:,1) = Ef;
end

% kxe(j) = kB*T/var(x(:,j));
% kye(j) = kB*T/var(y(:,j));


%% 
if strcmp(opts.output, 'tracks')
    tracks = cell(1,Nb);
    for j = 1:Nb
    	tracks{j} = [time',x(j,:)',y(j,:)', z(j,:)'];
    end

    varargout{1} = tracks;

    if nargout > 3
        msd = msdanalyzer(3, 'm','s','log');

        msd = msd.addAll(tracks);
        msd = msd.computeMSD;
        varargout{4} = msd;
    end

elseif strcmp(opts.output, 'data')
    data = struct('fName','BDsim','raw',struct(),'opts',opts);
    data.mPerPx = 1;
    data.raw.xCentresPx = x;
    data.raw.yCentresPx = y;
    data.raw.timeVecMs = time * 1e3;
    data = bead_preProcessCentres(data);
    if isfield(data.opts, 'Efreq')
        data.opts.Vfreq = data.opts.Efreq;
    end
    data.opts.gamma = gamma0;
    data.opts.D = D;
    varargout{1} = data;

    if nargout > 3
        data = bead_normMSD(data,'doPlots',false);
        varargout{4} = data.pro.amsdObj;
    end
end

if nargout > 1
    varargout{2} = struct('TotalForce',fs, 'ExternalForce', Efs);
end

if nargout > 2
    kxe = kB*T./var(x,0,2);
    kye = kB*T./var(y,0,2);
    kze = kB*T./var(z,0,2);
    varargout{3} = [kxe kye kze];
end


