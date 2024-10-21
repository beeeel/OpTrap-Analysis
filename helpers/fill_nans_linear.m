function filled_vector = fill_nans_linear(input_vector)
%% filled_vector = fill_nans_linear(input_vector)
% This function fills NaN values using linear interpolation between the
% non-NaN values on either side of each block of NaNs.

    % Create a copy of the input vector to modify
    filled_vector = input_vector;

    % Find the indices of NaNs
    nan_indices = isnan(input_vector);
    
    % Find the starting and ending indices of each block of NaNs
    diff_nan = diff([0; nan_indices; 0]);
    nan_starts = find(diff_nan == 1);
    nan_ends = find(diff_nan == -1) - 1;
    
    % Iterate over each block of NaNs
    for i = 1:length(nan_starts)
        start_idx = nan_starts(i);
        end_idx = nan_ends(i);
        
        % Determine the previous and next non-NaN values
        prev_idx = start_idx - 1;
        next_idx = end_idx + 1;
        
        % If there are no previous non-NaN values, use the next non-NaN value
        if prev_idx < 1
            filled_vector(start_idx:end_idx) = filled_vector(next_idx);
        % If there are no next non-NaN values, use the previous non-NaN value
        elseif next_idx > length(input_vector)
            filled_vector(start_idx:end_idx) = filled_vector(prev_idx);
        else
            % Perform linear interpolation between the surrounding non-NaN values
            filled_vector(start_idx:end_idx) = ...
                filled_vector(prev_idx) + ...
                (filled_vector(next_idx) - filled_vector(prev_idx)) * ...
                ((start_idx:end_idx) - prev_idx) / (next_idx - prev_idx);
        end
    end
end
