function opts = bead_loadOpts(data)

opts = struct();

% Load options file
if exist([data.dirPath '/opts.txt'],'file')
    fprintf('%s\t\t\tLoading options from file\nOptions Loaded:\n',repmat('#',30,1))
    fid = fopen([data.dirPath '/opts.txt'],'r');
    while ~feof(fid)
        ln = fgetl(fid);
        str = strsplit(ln,': ');
        try
            opts.(str{1}) = eval(str{2});
        catch ME
            warning('Caught: %s', ME.message)
        end
        fprintf('%s:\t\t\t%s\n',str{:})
    end
    fprintf('\n\n')
    fclose(fid);
end
