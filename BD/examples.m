%% How to use Brownian dynamics:
% The code is broken into two parts - one which creates an options
% structure that contains all the information needed to run a simulation,
% and the other which runs the simulation using the options in the
% structure. There are a number of default values which I've set up, and
% you can run simulations with different values as required. The simulation
% function can return data in two main ways: either as a cell of
% position-time tracks, or as a data structure suitable for using my
% microrheology code for further analysis and plotting.

%% Example 1: Single trapped bead

% 1) Create an opts structure that contains everything to do with the
% simulation. Parameters are given as name-value pairs, and the full list
% should be available using the command `help func_BDopts`.

eta = 0.87e-3; % Viscosity in Pa.s (defined here so it can be used later)
opts = func_BDopts('radius',3e-6,'kappaNm',[1 1 0.2] * 1e-6, ...
    'Nt', 5e5, 'eta', eta,'output','data');

% 2) Run the simulation and get the data struct out. This data struct is
% compatible with my other microrheology code (mostly found in
% OpTrap-Analysis/beads)

data = func_simulateDataField(opts);

% 3) Look at the data:

% Position-time data
bead_plotProData(data);
% Mean-squared displacement (a useful statistic to understand the
% position-time data)
data = bead_normMSD(data);

% 4) Fit to the position autocorrelation function to measure the viscosity
% and get trap stiffness using the variance of the data
data = bead_ACF(data);
% Trap relaxation time
tauc = data.pro.acfFit(1).fo.tauc;

Kb = 1.380649e-23;  % Boltzmann constant
T = 273 + 20;       % BD simulation defaults to temperature of 20 C
kappa = Kb * T ./ var(data.pro.xCentresM, 1, 2); % Trap stiffness

% Viscosity given by relaxation time, trap stiffness and bead size.
viscosity = tauc * kappa / 6 / pi / data.opts.radius;
fprintf('Calculated viscosity = %.4g Pa.s, expected %.4g Pa.s\n', viscosity, eta)

%% Example 2: Harmonically driven bead

% 0) Primer on anonymous functions. The harmonically driven case is built
% in if you use the input arguments Efreq and gamma, but here I'm going to
% define it explicitly to explain the use of anonymous function handles.

% Anonymous functions are functions defined like variables. First we need
% to define the constants used in the function:
E_amp = 1e3;    % Electric field amplitude in V/m
Efreq = 1;      % Driving frequency in Hz
% Then we declare the function by starting with an '@' and the input
% arguments in parenthesis. After that can be basically any matlab
% expression that returns an N x 3 x T matrix (N beads, T time points):
E_func = @(x, y, z, t) [E_amp 0 0] .* cos(2 * pi * Efreq * reshape(t, 1, 1, []));
% The '[E_amp 0 0]' means the force is in the X direction

% We can then use the anonymous function like a regular function
t = linspace(0,2,2e3); % Time array for input
E = E_func(0,0,0,t);   % The output will have size [1 x 3 x T] where T is the number of time points

figure
plot(t,squeeze(E(1,1,:)))
xlabel('Time (s)')
ylabel('E (V/m)')

% 1) Opts struct. This time we're simulating for 2 beads, set by the first
% input argument. We're providing an electric field and a valence for the
% bead (i.e. the number for q_bead is in units of proton charges)
opts = func_BDopts(2, 'radius',3e-6,'kappaNm',[1 1 0.2] * 1e-6, ...
    'E_func', E_func, 'q_bead', 1e3, ...
    'Nt', 5e5, 'eta', 0.87e-3, 'temp', 20, 'output','data');

% 2) Run the simulation. The multiple beads data all goes in one struct as
% different rows in the centres matrices (e.g. data.pro.xCentresM will have
% 2 rows and 500,000 columns)
data = func_simulateDataField(opts);

% 3) Look at the data:

% Position-time data
bead_plotProData(data);
% Mean-squared displacement (a useful statistic to understand the
% position-time data)
data = bead_normMSD(data);

% 4) The ACF can be analysed to extract the oscillation amplitude
data = bead_ACF(data);

% Oscillation ampltude in terms of the ratio of thermal to electrical forces
gamma = data.pro.acfFit(1).fo.gamma; 

%% Example 3: Double exponential relaxation

% 0) Define an anonymous function for the force. The double exponential 

% First we need to define the constants used in the function, two
% amplitudes and two rates, and one time offset (relaxation doesn't start
% at time 0).
tau = 1;    % Time (in seconds) to turn 'on' field
A0 = 1e4;   % Total amplitude (i.e. value at t=tau)
A1 = 7000;   % Amplitude of component 1
T1 = 0.25;  % Relaxation time of first component
T2 = 10;    % Relaxation time of second component
% Then we declare the function by starting with an '@' and the input
% arguments in parenthesis. 
E_func = @(x, y, z, t) [1, 1, 0] .* (reshape(t - tau,1,1,[]) >= 0) .* (A1 * exp(-reshape(t - tau,1,1,[]) ./ T1) ...
    + (A0 - A1) * exp(-reshape(t - tau,1,1,[]) ./ T2));
% The '[1, 1, 0]' means the force is in X and Y directions

% We can then use the anonymous function like a regular function
t = linspace(0,10,10e3); % Time array for input
E = E_func(0,0,0,t);   % The output will have size [1 x 3 x T] where T is the number of time points

figure
semilogy(t,squeeze(E(1,1,:)))
xlabel('Time (s)')
ylabel('E (V/m)')

% 1) Opts struct. This time we're simulating for 2 beads, set by the first
% input argument. We're providing an electric field and a valence for the
% bead (i.e. the number for q_bead is in units of proton charges)
opts = func_BDopts(1, 'radius',3e-6,'kappaNm',[1 1 0.2] * 1e-6, ...
    'E_func', E_func, 'q_bead', 1e3, ...
    'Nt', 5e5, 'eta', 0.87e-3, 'temp', 20, 'output','data');

% 2) Run the simulation. The multiple beads data all goes in one struct as
% different rows in the centres matrices (e.g. data.pro.xCentresM will have
% 2 rows and 500,000 columns)
data = func_simulateDataField(opts);

% 3) Look at the data:

% Position-time data
bead_plotProData(data);

% Now imagine you're going to do some more analysis elsewhere. To write the
% X position to a csv file, you can use this, but be careful because the
% file sizes can be very large due to the high sampling rate and low
% efficiency of the csv format. 
dlmwrite('example_X_position_data.csv', data.pro.xCentresM(1,:)','precision','%.9g');
% Note: Always use dlmwrite instead of csvwrite because csvwrite is limited
% to 5 significant figures, here I am writing 9 significant figures,
% defined by the precision argument.