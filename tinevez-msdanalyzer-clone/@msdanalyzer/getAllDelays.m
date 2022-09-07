function delays = getAllDelays(obj, indices)
% First, find all possible delays in time vectors.
% Time can be arbitrary spaced, with frames missings,
% non-uniform sampling, etc... so we have to do this clean.

if nargin < 2 || isempty(indices)
    indices = 1 : numel(obj.tracks);
end

n_tracks = numel(indices);
all_delays = cell(n_tracks, 1);
all_tSteps = cell(n_tracks, 1); 

% If times are duplicate, only do this once
if obj.time_duplicate
    n_tracks = 1;
end

for i = 1 : n_tracks
    index = indices(i);
    track = obj.tracks{index};
    t = track(:,1);
    tStep = msdanalyzer.roundn(mean(diff(t)), msdanalyzer.TOLERANCE);
    
    switch obj.sampling
        case 'linear'
            dT = tStep .* (0:size(t,1)-1)';
        case 'log'
            dT = tStep .* [0:10 ceil(obj.logbase.^(0:(size(t,1)-2)))];
            dT = dT( dT < max(t) )';
        case 'original'
            [T1, T2] = meshgrid(t, t);
            % Triu(A, 1) returns the lower half of A, not including the
            % main diagonal. Index to discard all 0 values.
            dT = msdanalyzer.roundn(abs(tril(T1 - T2, -1)), msdanalyzer.TOLERANCE);
            dT = dT(dT ~= 0);
            
            % % Original
            % dT = msdanalyzer.roundn(abs(T1(:)-T2(:)), msdanalyzer.TOLERANCE);
    end
    all_delays{i} = unique(dT);
    all_tSteps{i} = tStep;
end

if n_tracks > 1
    delays = unique( vertcat(all_delays{:}) );
else
    delays = all_delays{1};
end

% If tracks contains data not collected at the same time, all kinda havoc
% could happen!
if numel(unique([all_tSteps{:}])) ~= 1
    warning('Are time vectors equally spaced? Check yourself before you shrek yourself!')
end
end
