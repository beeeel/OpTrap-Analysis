function out = loadNyquist(filename)
% out = loadNyquist(filename)
%% Read csv from dropView and return nyquist data ([f, Z', Z''])
fid = fopen(filename);

str = '';

% Skip to the nyquist data
while ~contains(str, 'nyquist')
    str = fgetl(fid);
end
% Skip two more lines (blank line and table header)
[~] = fgetl(fid);
[~] = fgetl(fid);

out = NaN(3,1e3);
count = 1;
while ~feof(fid)
    str = fgetl(fid);
    try
        out(:,count) = sscanf(str, '"%g";"%g";"%g"');
        count = count + 1;
    catch
        break
    end
end
out = out([3 1 2], 1:count-1);