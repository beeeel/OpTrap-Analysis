function data = binaryToDouble(filename)
    % Load binary data from a file into a MATLAB double array.
    % Args:
    %     filename (string): Path to the binary file to be loaded.
    % Returns:
    %     data (double array): Array of doubles containing the data from the file.

    % Open the file in read-only mode with big-endian byte ordering
    fileID = fopen(filename, 'rb');
    
    if fileID == -1
        error('Failed to open file: %s', filename);
    end
    
    % Read the data as double precision floating-point numbers
    data = fread(fileID, Inf, 'double');
    
    % Close the file
    fclose(fileID);
end
