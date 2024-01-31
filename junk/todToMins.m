function mins = todToMins(tod, t0)
%% mins = todToMins(tod [, t0])
%% Convert time of day to minutes (optional: subtract initial time)
if ~exist('t0','var')
    t0 = 0;
end

if ~isnumeric(t0) || ~isscalar(t0) || ~isfinite(t0)
    error('Expected second input t0 to be finite numeric scalar')
end

hr0 = floor(t0/100);
mn0 = t0-100*hr0+60*hr0;

hrs = floor(tod/100);
mins = tod-100*hrs+60*hrs;

mins = mins - mn0;