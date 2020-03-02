function data = load_infos(treat, fieldList, varargin)
%%data = load_infos(treat, fieldList, [process])
%Load selected fields from saved info structs
%   Load every file in the infos folder that starts with 'info_$treat$'
% where $treat$ is the first input argument. From each file, take the
% data from every field in input argument fieldList, a cell vector.
%   Optional input argument - load differently processed datasets instead,
% e.g.: 'info-$process$_hela_$treat$' where $process$ is the third argin
%   Returns a (n+1 by m+1) cell array with a header column and filepath in
% first row, where n is the number of files found and m is the number of
% fields in fieldList

% Get the names of the files in the info directory
dir = '/home/ppxwh2/Documents/data/OpTrap/infos/';
fileList = strsplit(ls(dir));


% Count how many files we're gonna load - n_rows needs to be 1 larger
n_rows = 1;
for fileName = fileList
    % For each file name, split it by _, and if it's not empty then compare
    % the third element to the treat input ('info_hela_$treat$...')
    split = strsplit(fileName{:}, '_');
    % Also split it by - and if the process matches, ('info-$process$_hela_$treat$...')
    subsplit = strsplit(split{1},'-');
    if min(size([split{:}])) > 0 && strcmp(split{3}, treat) && (nargin == 3 && size(subsplit,2) > 1 && strcmp(subsplit{2}, varargin{1}))
        n_rows = n_rows + 1;
    end
end


% Preallocate a cell array - 1 row per file + 1 header, 1 column per field
% in fieldList + filepath field
data = cell(n_rows, 1 + size(fieldList, 2));
rowNum = 2;
data(1, :) = [{'filepath'}, fieldList];

% Actually load them and extract the fields
for fileName = fileList
    split = strsplit(fileName{:}, '_');
    subsplit = strsplit(split{1},'-');
    if min(size([split{:}])) > 0 && strcmp(split{3}, treat) && (nargin == 3 && size(subsplit,2) > 1 && strcmp(subsplit{2}, varargin{1}))
        % Load the file, place filepath in first column, and desired fields
        % in remaining columns
        disp(['Loading file: ', fileName])
        load(strcat(dir, fileName{:}),'info');
        data{rowNum, 1} = info(1).filepath;
        for idx = 1:max(size(fieldList,2))
            data{rowNum, idx+1} = [info.(fieldList{idx})];
        end
        rowNum = rowNum + 1;
    else
        disp(['Skipping file: ', fileName])
    end
end

end

