function [t, I, unit] = loadAmperometry(filename)
% [t, I, unit] = loadAmperometry(filename)
%% Load amperometry csv data from dropview

% filename = 'EIS/sample1_FA2_500mV_precond-500mV_10s.csv';

s = importdata(filename);
str = strsplit(s{3}, {'Current (', ')'});
unit = str{3};
s = s(4:end);

t = nan(size(s));
I = nan(size(s));
for ind = 1:length(s)
    str         = strsplit(s{ind},{'";"','"'});
    t(ind)      = str2double(str{2});
    I(ind)      = str2double(str{3});
end