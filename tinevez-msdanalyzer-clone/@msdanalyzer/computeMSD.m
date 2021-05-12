function obj = computeMSD(obj, indices)
%%COMPUTEMSD Compute the mean-squared-displacement for this object.
%
% obj = obj.computeMSD computes the MSD for all the tracks stored
% in this object. If a drift correction was computed prior to this
% method call, it is used to correct positions before MSD
% calculation.
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

n_tracks = numel(indices);
fprintf('Computing MSD of %d tracks... ', n_tracks);

% First, find all possible delays in time vectors.
% Time can be arbitrary spaced, with frames missings,
% non-uniform sampling, etc... so we have to do this clean.
% We use a certain tolerance to bin delays together
delays = obj.getAllDelays;
n_delays = numel(delays);

fprintf('Computing MSD using %d delays... ', n_delays);

obj.msd = cell(n_tracks, 1);
if ~isempty(obj.drift)
    tdrift = obj.drift(:,1);
    xdrift = obj.drift(:, 2:end);
end

fprintf('%5d/%5d', 0, n_tracks);

for i = 1 : n_tracks
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d/%5d', i, n_tracks);
      
    mean_msd    = zeros(n_delays, 1);
    M2_msd2     = zeros(n_delays, 1);
    n_msd       = zeros(n_delays, 1);
    
    index = indices(i);
    track = obj.tracks{index};
    t = track(:,1);
    t = msdanalyzer.roundn(t, msdanalyzer.TOLERANCE);
    X = track(:, 2:end);
    
    n_detections = size(X, 1);

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
    
    
    
    if sum(strcmp(obj.sampling, {'linear', 'original'})) == 1
            for j = 1 : n_detections - 1
                
                switch obj.sampling
                    case 'linear'
                        index_in_all_delays = 2 : n_detections - j + 1;
                    case 'log'
                        % Don't think this will ever run
                        warning('I didn''t think this will ever run')
                        n = 0 : n_detection - j - 1;
                        index_in_all_delays = 1 + ceil(obj.logbase.^n);
                    case 'original'
                        % Delay in physical units
                        dt = t(j+1:end) - t(j);
                        dt = msdanalyzer.roundn(dt, msdanalyzer.TOLERANCE);
                        
                        % Determine target delay index in bulk
                        n_dt = numel( dt );
                        index_in_all_delays = NaN( n_dt, 1 );
                        
                        for k = 1 : n_dt
                            [ ~, l ] = min( abs( dt(k) - delays ) );
                            index_in_all_delays( k ) = l;
                        end
                end
                
                % Square displacement in bulk
                dX = X(j+1:end,:) - repmat(X(j,:), [(n_detections-j) 1] );
                dr2 = sum( dX .* dX, 2); % Each dimension has 1 column
                
                % Store for mean computation / Knuth
                n_msd(index_in_all_delays)     = n_msd(index_in_all_delays) + 1;
                delta = dr2 - mean_msd(index_in_all_delays);
                mean_msd(index_in_all_delays) = mean_msd(index_in_all_delays) + delta ./ n_msd(index_in_all_delays);
                M2_msd2(index_in_all_delays)  = M2_msd2(index_in_all_delays) + delta .* (dr2 - mean_msd(index_in_all_delays));
            end
    elseif strcmp(obj.sampling, 'log')
            %             n = 0 : n_detection - 1;
            %             delays = ceil(obj.logbase.^n);
            alldTs = unique([0 ceil(obj.logbase.^(0:(size(t,1)-2)))]);
            alldTs = alldTs( alldTs < size(t,1) );
            for j = 2:length(alldTs)
                dT = alldTs(j);
                % Calculate all square displacements for this delay
                dX = X(dT+1:end,:) - X(1:end-dT,:);
                dr2 = sum( dX .* dX, 2); % Each dimension has 1 column
                
                % Store
                n_msd(j) = size(dr2,1);
                mean_msd(j) = mean(dr2,1);
                M2_msd2(j) = var(dr2,1,1)*n_msd(j);
            end 
    end
    n_msd(1) = n_detections;
    std_msd = sqrt( M2_msd2 ./ n_msd ) ;
    
    % We replace points for which N=0 by Nan, to later treat
    % then as missing data. Indeed, for each msd cell, all the
    % delays are present. But some tracks might not have all
    % delays
    delay_not_present = n_msd == 0;
    mean_msd( delay_not_present ) = NaN;
    
    obj.msd{index} = [ delays mean_msd std_msd n_msd ];
    
end
fprintf('\b\b\b\b\b\b\b\b\bDone.\n')

obj.msd_valid = true;

end