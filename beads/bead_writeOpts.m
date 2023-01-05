function bead_writeOpts(data, opts)

optsName = [data.dirPath '/opts.txt'];

if ~exist(optsName,'file') || getgood(sprintf('Overwrite %s?',optsName))
    fid2 = fopen(optsName,'w');
    fnames = fieldnames(opts)';
    for fld = fnames
        val = eval(sprintf('opts.%s',fld{1}));
        switch fld{1}
            case 'skipSuffixes'
                fprintf(fid2,'%s: %s%s%s\n',fld{1},'''',val,'''');
            otherwise
                fprintf(fid2,'%s: %s\n',fld{1}, val);
%                 if isscalar(val)
%                     fprintf(fid2,'%s: %i\n',fld{1},eval(sprintf('opts.%s',fld{1})));
%                 else
%                     fprintf(fid2,'%s: [', fld{1});
%                     for idx = 1:size(val,1)
%                         for jdx = 1:size(val,2)
%                             fprintf(fid2,' %i ',val(idx,jdx));
%                         end
%                         fprintf(fid2,';');
%                     end
%                     fprintf(fid2,']\n');
%                 end
        end
    end
    fclose(fid2);
end

end


function result = getgood(message)
failure = true;
while failure
    try
        result = logical(input(message));
        failure = false;
    catch ME
        warning(ME.message)
        fprintf('\nTry again! 0 or 1\n')
    end
end
end