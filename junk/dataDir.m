function dir = dataDir(subDir)
%% Given a subdirectory, return a machine-specific full path for data
[~, hn]  = system('hostname');
if strncmp(hn, 'will-linux',10)
    dir = ['/home/will/Data/work/' subDir];
else
    dir = ['/home/scan/will/' subDir];
end