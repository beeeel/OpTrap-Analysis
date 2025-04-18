function opts = func_BDopts(varargin)
%% Build an opts struct for BD sims
% opts = func_BDopts([Nbeads], varargin)
% Supply options in name-value pairs, any omitted will be set to defaults.
% Possible options, defaults, and units:
%   radius      2.5e-6          m
%   temp        20              °C
%   pos0        [0 0 0]         m
%   kappaNm     [1 1 .3]*1e-6   N/m
%   E_func      []              V/m  (old default: @(x,y,z,t)0 )
%   F_func      []              N  (old default: @(x,y,z,t)0 )
%   q_bead      1               e (fundamental charge, 1.6e-19 C)
%   dt          1e-4            s
%   Nt          1e5             samples
%   eta         0.97e-3         Pas
%   rng_seed    [               dimensionless
%   output      'data'          'tracks' OR 'data'
%   m_bead      0               kg
%
% Particle tracking errors can be simulated
%   statErrM2   []              m^2 (Variance of apparent position for stationary bead)
%   dynErrS     []              s (Effective exposure time, must be much larger than dt)
%
% Instead of trap stiffness, you may choose to give corner time
%   tauc        []              s
%
% Instead of defining a function, you may choose a value of gamma (force
% ratio) and frequency - this will define an E field in X only
%   gamma       []              dimensionless
%   Efreq       []              Hz

if isnumeric(varargin{1})
    N_beads = varargin{1};
    args = varargin(2:end);
else
    N_beads = 1;
    args = varargin;
end

posNumeric = @(x) isnumeric(x) && (size(x,1) == N_beads || size(x,1) == 1 ) && all(x > 0);
Numeric = @(x) isnumeric(x) && (size(x,1) == N_beads || size(x,1) == 1 );
posVec =  @(x) isnumeric(x) && (size(x,1) == N_beads || size(x,1) == 1 ) && size(x,2) == 3;
posScalar = @(x) isnumeric(x) && isscalar(x) && x > 0;

p = inputParser;

p.addOptional('radius',  2.5e-6,        posNumeric);
p.addOptional('temp',    20,            Numeric);
p.addOptional('pos0',    [0 0 0],       posVec);
p.addOptional('kappaNm', [1 1 .3]*1e-6, posVec);
p.addOptional('E_func',  [],  @(x) isa(x, 'function_handle') && isnumeric(x(0,0,0,0))); % old default @(x,y,z,t) zeros(size(x,1),3,numel(t))
p.addOptional('F_func',  [],  @(x) isa(x, 'function_handle') && isnumeric(x(0,0,0,0))); % old default @(x,y,z,t) zeros(size(x,1),3,numel(t))
p.addOptional('q_bead',  1,             posNumeric);
p.addOptional('m_bead',  0,             posNumeric);
p.addOptional('dt',      1e-4,          posScalar);
p.addOptional('Nt',      1e5,           posScalar);
p.addOptional('eta',     0.97e-3,       posNumeric); % % eta = 2.414E-5 * 10.^(247.8 ./(T -140)); 
p.addOptional('rng_seed',[],             @(x) isscalar(x) && isinteger(x) && x < 2^32);
p.addOptional('output',  'data',        @(x) any(strcmp(x, {'tracks','data'})));

p.addOptional('statErrM2',[],           posScalar);
p.addOptional('dynErrS',  [],           posScalar);

p.addOptional('tauc',    [],            posNumeric);

p.addOptional('Efreq',   [],            posScalar);
p.addOptional('gamma',   [],            posScalar);

parse(p, args{:});

opts = p.Results;

if ~isempty(opts.tauc)
    if posScalar(opts.tauc)
        opts.kappaNm = [1 1 .3] * (6 * pi * opts.radius * opts.eta) ./ opts.tauc;
    else
        opts.kappaNm = (6 * pi * opts.radius * opts.eta) ./ opts.tauc;
    end
else
    opts.tauc = (6 * pi * opts.radius * opts.eta) ./ opts.kappaNm;
end

if ~isempty(opts.Efreq) && ~isempty(opts.gamma)
    fc = opts.kappaNm(1) ./ (6*pi*opts.radius*opts.eta);
    E_amp =  sqrt(opts.gamma * 2 * opts.kappaNm(1)*kBT(opts.temp+273) * (1 + opts.Efreq^2/fc^2)/(opts.q_bead*1.6e-19)^2);          % V/m

    opts.E_func = @(x, y, z, t) [E_amp 0 0] .* cos(2 * pi * opts.Efreq * reshape(t, 1, 1, []));
end

fns = fieldnames(opts);

for idx = 1:length(fns)
    if ~strcmp(fns{idx}, {'E_func','F_func', 'Nt','dt','rng_seed', 'output','Efreq','gamma'})
        if size(opts.(fns{idx}), 1) == 1
            opts.(fns{idx}) = repmat(opts.(fns{idx}), N_beads, 1);
        elseif size(opts.(fns{idx}),1) ~= N_beads && ~isempty(opts.(fns{idx}))
            error('Unexpected 1st dim size for options %s: %i', fns{idx}, size(opts.(fns{idx}),1))
        end
    end
end



opts.Nbeads = N_beads;

end
