function obj = computeDdist(obj, indices, normR, dTs)
%%COMPUTEDdist Compute the displacement distribution for this object.
%
% obj = obj.computeDdist computes the displacement distribution for all the
% tracks stored in this object. If a drift correction was computed prior to
% this method call, it is used to correct positions before MSD calculation.
%
% obj = obj.computeDdist(indices) computes the Ddist only for the
% particles with the specified indices. Use an empty array to take
% all particles.
%
% obj = obj.computDdist([], normD) computes Ddist divided by normalisation
% factor normD. Be careful to match the units of normD to the spatial units
% in the msdanalyzer!
%
% obj = obj.computDdist([], [], dTs) computes Ddist at a user-specified
% set of lag times, if they are available in the data. Be careful to match
% the units of dTs to the temporal units in the msdanalyzer!
% 
% Output: obj.Ddist, a cell array with one row per track. First element
% contains the lag times queried, the number of observations, and the NGP.
% Second element contains histogram counts and bin edges, such that columns
% 1:floor(end/2) are the counts and ceil(end/2):end are the edges.
%
% NGP - Bursac 2005, Weeks 2002, Rahman 1964.
%  <z^4> / (3 <z^2> ^2 ) - 1
ngp = @(z, rho) sum(rho.*(z.^4),1) ./ (3 .* sum(rho.*(z.^2),1) .^2) - 1;

if nargin < 2 || isempty(indices)
    indices = 1 : numel(obj.tracks);
end

if obj.n_dim > 1
    error('Ddist only works in 1D')
end

if exist('normR','var') && ~isempty(normR) && numel(indices) ~= numel(normR)
    error('Need 1 normR for each track')
end



n_tracks = numel(indices);
fprintf('Computing Ddist of %d tracks... ', n_tracks);

%fprintf('Computing Ddist using %d delays... ', n_delays);

%obj.msd = cell(n_tracks, 1);
if ~isempty(obj.drift)
    tdrift = obj.drift(:,1);
    xdrift = obj.drift(:, 2:end);
end

fprintf('%5d/%5d', 0, n_tracks);
warning('off','stats:lillietest:OutOfRangePLow')
for i = 1 : n_tracks
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d/%5d', i, n_tracks);
      
    index = indices(i);
    track = obj.tracks{index};
    t = track(:,1);
    t = msdanalyzer.roundn(t, msdanalyzer.TOLERANCE);
    X = track(:, 2:end);
        
    dt = t(2) - t(1);
    
    %%% Calculate indices (alldTs needs to be in indexing units!)
    
    % Don't need this if only doing preset ΔD    lb = 2; % Logbase
    % A bit overkill!    alldTs = unique([0 ceil(lb.^(0:(size(t,1)-2)))])/dt;
    
    minInd = 10; % Minimum independent observations of maxDelay
    
    % Take the dTs from the object unless supplied by user.
    if ~exist('dTs', 'var') || isempty(dTs)
        dTs = obj.dTs;
    end
    alldTs = round(dTs / dt); % Ends up in indexing units!
    alldTs = alldTs(alldTs ~= 0);
    alldTs = alldTs( alldTs < size(t,1)/minInd )';
    
    % Histogram bins - centres then edges
    z = -7:0.1:7;
    bw = (z(2)-z(1))/2;
    edg = [z(1)-bw, bw + z]';
%     edg = 0.05 + (-20:0.1:20)';
    nBins = numel(edg) - 1;
    
    if ~exist('normR', 'var') || isempty(normR)
        norm = @(x) (x - mean(x,'all')) ./ std(x);
    else
        norm = @(x) (x - mean(x,'all')) ./ normR(i);
    end
    
    n_delays    = numel(alldTs);
    counts      = zeros(nBins, n_delays);
    edges       = repmat(edg, 1, n_delays);
    n_msd       = zeros(n_delays, 1);
%     pvals       = zeros(n_delays, 2);
    alpha       = zeros(n_delays, 1);

    % Determine drift correction
    if ~isempty(obj.drift)
        % Determine target delay index in bulk
        [~, index_in_drift_time, index_in_track_time] = intersect(tdrift, t);
        % Keep only track times that can be corrected.
        X = X(index_in_track_time, :);
        % Subtract drift position to track position
        X = X - xdrift(index_in_drift_time, :);
        
    end
    
    % Calculate increment distribution
    for j = 1:n_delays
        dT = alldTs(j);
        % Calculate all displacements for this delay
        dX = X(dT+1:end,:) - X(1:end-dT,:);
        dX = norm(dX); % Perform normalisation as per Bursac et al 2005
        
        % Bin like a histogram
        [N, ~] = histcounts(dX, edg);%, 'Normalization', 'probability');
        
        n_msd(j) = size(dX,1);
        counts(:,j) = N;
        
        alpha(j) = ngp(z', N'./sum(N));
%         [~, pvals(j,1)] = ttest(dX);
%         [~, pvals(j,2)] = lillietest(dX);
    end
    
    
    obj.Ddist{index,1} = [ alldTs*dt n_msd alpha];% pvals ];
    % First cell contains the increment times queries, number of points ngp
    
    obj.Ddist{index,2} = [ counts; edges ]; % Sorry about how this looks.
    % This way to access the counts and edges for a given dT, you get the
    % fast slice indexing and predictable sizes: 2×nBins+1 by n_delays
    
end
fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
warning('on','stats:lillietest:OutOfRangePLow')

end