function addpath_recursive(folder)

ignoreDirs = {'.','..','.git'};
addpath(folder)

subdirs = dir(folder);
for obj = subdirs'
    if obj.isdir && ~any(strcmp( obj.name, ignoreDirs)) && ~strcmp(obj.name(1),'@')
        addpath_recursive(strcat(folder, obj.name, '/'))
    end
end

end
