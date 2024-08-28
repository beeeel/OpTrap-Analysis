function filled_vector = fill_nans_last_value(input_vector)
%% filled_vector = fill_nans_last_value(input_vector)
% an idea for phase unwrapping in signal reconstruction

    % Create a copy of the input vector to modify
    filled_vector = input_vector;
    
    % Iterate through the vector
    for i = 2:length(input_vector)
        if isnan(filled_vector(i))
            % Replace NaN with the last non-NaN value
            filled_vector(i) = filled_vector(i-1);
        end
    end
end