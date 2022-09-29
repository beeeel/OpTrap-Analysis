function obj = computeDdist(obj, indices, normD)
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

if nargin < 2 || isempty(indices)
    indices = 1 : numel(obj.tracks);
end

if obj.n_dim > 1
    error('Ddist only works in 1D')
end

if ~exist('normD','var')
    normD = ones(size(indices));
end

if numel(indices) ~= numel(normD)
    error('Need 1 normD for each track')
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
    
    minInd = 1; % Minimum independent observations of maxDelay
    
    % Take the dTs from the object.
    alldTs = round(obj.dTs / dt); % Ends up in indexing units!
    alldTs = alldTs( alldTs < size(t,1)/minInd )';
    
    % Number of histogram bins
    edg = 0.05 + (-7:0.1:7)';
    nBins = numel(edg) - 1;
    
    n_delays    = numel(alldTs);
    counts      = zeros(nBins, n_delays);
    edges       = repmat(edg, 1, n_delays);
    n_msd       = zeros(n_delays, 1);
%     pvals       = zeros(n_delays, 2);
    %covs        = cell(n_delays, 1);

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
        norm = @(x) (x - mean(x,'all')) ./ std(x);
        dX = norm(dX); % Perform normalisation as per Bursac et al 2005
        
        % Bin like a histogram
        [N, ~] = histcounts(dX, edg);%, 'Normalization', 'probability');
        
        n_msd(j) = size(dX,1);
        counts(:,j) = N;
%         [~, pvals(j,1)] = ttest(dX);
%         [~, pvals(j,2)] = lillietest(dX);
    end
    
    
    obj.Ddist{index,1} = [ alldTs*dt n_msd];% pvals ];
    % First cell contains the increment times queries, number of points and
    % probabilities that 1) Mean == 0, 2) Distribution == normal.
    
    obj.Ddist{index,2} = [ counts; edges ]; % Sorry about how this looks.
    % This way to access the counts and edges for a given dT, you get the
    % fast slice indexing and predictable sizes: 2×nBins+1 by n_delays
    
end
fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
warning('on','stats:lillietest:OutOfRangePLow')

end