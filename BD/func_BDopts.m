function opts = func_BDopts(varargin)
%% Build an opts struct for BD sims
% opts = func_BDopts([Nbeads], varargin)
% Supply options in name-value pairs, any omitted will be set to defaults.
% Possible options, defaults, and units:
%   radius      2.5e-6          m
%   temp        20              °C
%   pos0        [0 0 0]         m
%   kappaNm     [1 1 .3]*1e-6   N/m
%   E_func      @(x,y,z,t)0     V/m
%   q_bead      1               e (fundamental charge, 1.6e-19 C)
%   dt          1e-4            s
%   Nt          1e5             samples
%   eta         0.97e-3         Pas
%   rng_seed    0               dimensionless

if isnumeric(varargin{1})
    N_beads = varargin{1};
    args = varargin(2:end);
else
    N_beads = 1;
    args = varargin;
end

posNumeric = @(x) isnumeric(x) && (size(x,1) == N_beads || size(x,1) == 1 ) && all(x > 0);
posVec =  @(x) isnumeric(x) && (size(x,1) == N_beads || size(x,1) == 1 ) && size(x,2) == 3;
posScalar = @(x) isnumeric(x) && isscalar(x) && x > 0;

p = inputParser;

p.addOptional('radius',  2.5e-6,        posNumeric);
p.addOptional('temp',    20,            posNumeric);
p.addOptional('pos0',    [0 0 0],       posVec);
p.addOptional('kappaNm', [1 1 .3]*1e-6, posVec);
p.addOptional('E_func',  @(x,y,z,t) 0,  @(x) isa(x, 'function_handle') && isnumeric(x(0,0,0,0)));
p.addOptional('q_bead',  1,             posNumeric);
p.addOptional('dt',      1e-4,          posScalar);
p.addOptional('Nt',      1e5,           posScalar);
p.addOptional('eta',     0.97e-3,       posNumeric);
p.addOptional('rng_seed',1,             @(x) isscalar(x) && isinteger(x) && x < 2^32);

parse(p, args{:});

opts = p.Results;

fns = fieldnames(opts);

for idx = 1:length(fns)
    if ~strcmp(fns{idx}, {'E_func', 'Nt','dt','rng_seed'})
        if size(opts.(fns{idx}), 1) == 1
            opts.(fns{idx}) = repmat(opts.(fns{idx}), N_beads, 1);
        elseif size(opts.(fns{idx}),1) ~= N_beads
            error('Unexpected 1st dim size for options %s: %i', fns{idx}, size(opts.(fns{idx}),1))
        end
    end
end

opts.Nbeads = N_beads;

end