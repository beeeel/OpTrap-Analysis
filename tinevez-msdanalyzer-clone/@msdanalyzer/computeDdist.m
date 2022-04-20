function obj = computeDdist(obj, indices)
%%COMPUTEDdist Compute the displacement distribution for this object.
%
% obj = obj.computeDdist computes the displacement distribution for all the
% tracks stored in this object. If a drift correction was computed prior to
% this method call, it is used to correct positions before MSD calculation.
%
% Results are stored in the msd field of this object as a cell
% array, one cell per particle. The array is a double array of size
% N x 4, and is arranged as follow: [dt mean std N ; ...] where dt
% is the delay for the MSD, mean is the mean MSD value for this
% delay, std the standard deviation and N the number of points in
% the average.
%
% obj = obj.computeMSD(indices) computes the MSD only for the
% particles with the specified indices. Use an empty array to take
% all particles.

if nargin < 2 || isempty(indices)
    indices = 1 : numel(obj.tracks);
end

if obj.n_dim > 1
    error('Ddist only works in 1D')
end

n_tracks = numel(indices);
fprintf('Computing Ddist of %d tracks... ', n_tracks);

%fprintf('Computing Ddist using %d delays... ', n_delays);

obj.msd = cell(n_tracks, 1);
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
    
    n_detections = size(X, 1);
    
    % Calculate indices (alldTs needs to be in indexing units!)
    lb = 2; % Logbase
    minInd = 10; % Minimum independent observations of maxDelay
    
%     alldTs = unique([0 ceil(lb.^(0:(size(t,1)-2)))])/dt;
    % Actually, just do these 4 (if available).
    alldTs = round([1e-3 1e-2 1 10 50]/dt);
    alldTs = alldTs( alldTs < size(t,1)/minInd )';
    
    % Number of histogram bins
    nBins = 100;
    
    n_delays    = numel(alldTs);
    counts      = zeros(nBins, n_delays);
    edges       = zeros(nBins + 1, n_delays);
    n_msd       = zeros(n_delays, 1);
    pvals       = zeros(n_delays, 2);
    covs        = cell(n_delays, 1);

    % Determine drift correction
    if ~isempty(obj.drift)
        % Determine target delay index in bulk
        [~, index_in_drift_time, index_in_track_time] = intersect(tdrift, t);
        % Keep only track times that can be corrected.
        X = X(index_in_track_time, :);
        t = t(index_in_track_time);
        % Subtract drift position to track position
        X = X - xdrift(index_in_drift_time, :);
        
    end
    
    % Calculate increment distribution
    for j = 1:n_delays
        dT = alldTs(j);
        % Calculate all square displacements for this delay
        dX = X(dT+1:end,:) - X(1:end-dT,:);
        
        % Bin like a histogram
        [N, edg] = histcounts(dX, nBins);
        
        % Calculate autocorrelation (covariance)
        C = xcorr(dX);
        tau = dt*(1:ceil(size(C,1)/2))'; % Delay time for covariance
        
        n_msd(j) = size(dX,1);
        counts(:,j) = N;
        edges(:,j) = edg;
        [~, pvals(j,1)] = ttest(dX);
        [~, pvals(j,2)] = lillietest(dX);
        covs{j} = [tau, C(ceil(end/2):end)];
    end
    
    
    obj.Ddist{index,1} = [ alldTs*dt n_msd pvals ];
    % First cell contains the increment times queries, number of points and
    % probabilities that 1) Mean == 0, 2) Distribution == normal.
    obj.Ddist{index,2} = [ counts; edges ]; % Sorry about how this looks.
    % This way to access the counts and edges for a given dT, you get the
    % fast slice indexing and predictable sizes: 2×nBins+1 by n_delays
    obj.Ddist{index,3} = covs;
    % Autocorrelation vectors are in cells because different sizes.
end
fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
warning('on','stats:lillietest:OutOfRangePLow')
obj.msd_valid = true;

end