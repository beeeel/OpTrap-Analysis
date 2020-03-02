function [ Imstack ] = load_imstack( dir, treat, speed, trap, spots, rep)
%Imstack = load_imstack( dir, treat, speed, trap, spots, rep) - wrapper for bfopen to open OME/TIFF files
%   Takes inputs as strings. Use numeric 0 to skip an argument
%   Input arguments:
%       dir     - directory relative to OpTrap where data is
%       treat   - treatment cell received, e.g. ctrl
%       speed   - speed in 100 * mm/s (filename parameter)
%       trap    - trap strength as % max laser power
%       spots   - number of spots used. Set as 0 to leave out of filename
%       rep     - repeat at those parameters. Set as 0 to leave out.
%   Output argument:
%       Imstack - cell array loaded from OME/TIFF file
%
%   Example filename: OpTrap/0508/hela_ctrl_s_040_tr_70_MMStack.ome.tif
%       dir = '0508',     treat = 'ctrl',   speed = '040',    trap = '70', 
%       spots = 0, rep = 0.

% Construct the file path variable
path = strcat('/home/ppxwh2/Documents/data/OpTrap/');
if dir ~= 0
    path = strcat(path, dir, '/');
end

addpath(path);

file = strcat('hela_', treat, '_s_', speed, '_tr_', trap);


if spots ~= 0
    file = strcat(file, '_sp_', spots);
end
if rep ~= 0
    file = strcat(file, '_', rep);
end

filepath = strcat(path, file, '_MMStack.ome.tif');
disp(strcat('Checking file: ',filepath));

if exist(filepath, 'file') == 2
    % bfopen(filename) - the bioformats file loader. Returns a cell
    Imstack = bfopen(filepath);
else
    % All this work for a nice error message
    disp('Did you remember to move the stacks out of their sub-folders?');
    error(['File does not exist - check the directory/filename is '...
        'correct and try again? Filenames are case sensitive']);
end


end

