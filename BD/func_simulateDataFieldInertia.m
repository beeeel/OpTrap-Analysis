function [varargout] = func_simulateDataFieldInertia(opts)
%% [tracks, Xforces, stiffnesses, msdanalyzer] = func_simulateDataFieldInertia(opts)
%% BD of particle in optical trap and electric field, inertial regime
% Input an opts structure from func_BDopts, outputs optional - tracks for
% msdanalyzer, stiffnesses (kBT/variance) as [kx1 ky1 kz1; kx2...],
% msdanalyzer with msd calculated.
%
% Based on Volpe 2013, adapted to include external force by Will Hardiman,
% August 2024

Nb = opts.Nbeads;

r = opts.radius;                       % radius of bead (m)
m = opts.m_bead;                         % Mass of bead (kg)
t = opts.temp;                         % temperature in Celsius
T = 273+t;                      % temperature in Kelvin

kx = opts.kappaNm(:,1);                      % trap stiffness along x
ky = opts.kappaNm(:,2);                      % trap stiffness along y
kz = opts.kappaNm(:,3);                      % trap stiffness along y

if isfield(opts, 'E_func') && ~isempty(opts.E_func)
    f = opts.E_func;		% Force field defined via function of position
    fSF = 1.6e-19 * opts.q_bead;			% Charge of particle in Coulombs
elseif isfield(opts, 'F_func') && ~isempty(opts.F_func)
    f = opts.F_func;
    fSF = 1;
else
    f = @(x,y,z,t) zeros(size(x,1),3,length(t));
    fSF = 1;
end

eta = opts.eta;			% Viscosity of solution

%viscosity-temperature eq for water 
%eta = 2.414E-5 * 10.^(247.8 ./(T -140)); 
%complex viscosity = dynamic viscosity + out-of-phase viscosity (for
%water = dynamic viscosity)

N = opts.Nt;                       % number of samples
deltat = opts.dt;                 % duration of each frame
fps = 1/deltat;                      % sample rate
gamma0 = 6*pi*r.*eta;            % viscous drag force


kB=1.38E-23;                    % Boltzmann's constant
D=(kB.*T)./gamma0;                  % diffusion coefficient

% Check timestep vs trap stiffness
if any(gamma0./opts.kappaNm < 10/fps)
    warning('Large time step compared to trap relaxation time')
end
% initialise bead position vectors
pos0 = opts.pos0;
x = zeros(Nb,N)+pos0(:,1);
y = zeros(Nb,N)+pos0(:,2);
z = zeros(Nb,N)+pos0(:,3);
fxs = zeros(Nb, 3, N);
fEs = zeros(Nb, 3, N);
time = (0:N-1).*opts.dt;

% random Gaussian noise
if opts.rng_seed ~= 0; rng(opts.rng_seed); end
noise = randn(Nb, 3, N) .* sqrt(2 * kB * T .* gamma0 / deltat) ;

%% Simulation

for i=2:N-1
    
    
    fE = f(x(:,i), y(:,i), z(:,i), time(i)) .* fSF;
    
    % Calculate position (Langevin's equation) - Volpe and Volpe 2012,
    % combining elements of eq 6 and 15, then solving as in eq 7+8,
    % including external force (my contribution)
    x(:,i+1) = 1 ./ (m / deltat^2 + gamma0 / deltat) .*      ... % Inertial and drag scale factor
        (x(:,i) .* (2 .* m / deltat^2 + gamma0/deltat)       ... % Inertial, and drag force
        - (x(:,i) - pos0(1)) .* kx                           ... % Trap force
        + x(:,i-1) .* (-m / deltat^2) + noise(:,1,i) + fE(:,1)); % Inertia, noise, and external force
    
    y(:,i+1) = 1 ./ (m / deltat^2 + gamma0 / deltat) .*      ... % Inertial + drag scale factor
        (y(:,i) .* (2 .* m / deltat^2 + gamma0/deltat)       ... % Inertial, and drag force
        - (y(:,i) - pos0(2)) .* ky                           ... % Trap force
        + y(:,i-1) .* (-m / deltat^2) + noise(:,2,i) + fE(:,2)); % Inertia, noise, and external force
    
    z(:,i+1) = 1 ./ (m / deltat^2 + gamma0 / deltat)         ... % Inertial + drag scale factor
        .* (z(:,i) .* (2 .* m / deltat^2 + gamma0/deltat)    ... % Inertial, and drag force
        - (z(:,i) - pos0(3)) .* kz                           ... % Trap force
        + z(:,i-1) .* (-m / deltat^2) + noise(:,3,i) + fE(:,3)); % Inertia, noise, and external force
    
    fxs(:,:,i) = [- (x(:,i) - pos0(:,1)) .* kx, ...	% Trap force
        noise(:,1,i),                      ...  % Diffusion
        fE(:,1)];                               % External force
    
    fEs(:,:,i) = fE;
end

% Simulate dynamic error by downsampling position and time traces.
addDynamicError();

% Simulate static error by adding Gaussian noise to position data (not
% force or time)
addDynamicError();

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
    varargout{2} = struct('TotalForce',fxs, 'ExternalForce', fEs);
end

if nargout > 2
    kxe = kB*T./var(x,0,2);
    kye = kB*T./var(y,0,2);
    kze = kB*T./var(z,0,2);
    varargout{3} = [kxe kye kze];
end

    function addDynamicError()
        if isfield(opts, 'dynErrS') && ~isempty(opts.dynErrS) && opts.dynErrS > 0
            N = opts.dynErrS / opts.dt;
            if N ~= round(N)
                N = round(N);
                warning('Changed exposure time (dynamic error) to %g s', N*opts.dt)
            end
            order1 = [2 3 1];
            order2 = [3 2 1]; %#ok<NASGU>
            for fn = ["x", "y","z","time", "fs", "Efs"]
                tmp = permute(eval(fn), order1);
                tmp = reshape(tmp, N, [], size(tmp,3));
                tmp = mean(tmp, 1); %#ok<NASGU>
                eval(sprintf('%s = ipermute(tmp, order2);', fn))
            end
        end
    end

    function addStaticError()
        if isfield(opts, 'statErrM2') && ~isempty(opts.statErrM2) && opts.statErrM2 > 0
            for fn = ["x", "y", "z"]
                tmp = randn(Nb,N); %#ok<NASGU>
                eval(sprintf('%s = %s + tmp;', fn, fn))
            end
        end
    end

end
