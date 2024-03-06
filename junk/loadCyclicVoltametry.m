function [V, I, units] = loadCyclicVoltametry(filename)
% [t, I, unit] = loadAmperometry(filename)
%% Load amperometry csv data from dropview

% filename = '../2024_01_26/NHcell_CV2_pm1500mV_100mVs_3scans.csv';

s = importdata(filename);
str = strsplit(s{3}, {'Potential (', ')', 'Current (', ')'});
units = str([2 4]);
s = s(4:end);

V = nan(size(s));
I = nan(size(s));
for ind = 1:length(s)
    try
        str         = strsplit(s{ind},{'";"','"'});
        if ~isempty(str{2})
            V(ind)      = str2double(str{2});
            I(ind)      = str2double(str{3});
        end
    catch ME
        error('Loading probably failed because the files are written as multiple experiments each consisting 1 scan. Fix your code!')
        break
        error(ME.message)
    end
end