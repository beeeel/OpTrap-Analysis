%% Cut up a dat file to make a bunch of new ones
% Because I was foolish when I programmed the MSD code and it can't handle
% an array of arbritrary cropTs.

dataPath = '../2022_11_03/cond1_stretcher1/disc0_bead_5/rep1_3_laser_off';

cropT = [1          1918    % These values were found by much pain and 
         1920       3966    % will be passed down from generation to 
         3968       6014    % generation to honour the difficulty I went 
         6016       8060    % through trying to find how to get data from 
         8065       10000]; % a data cursor on a plot only to find that I 
                            % couldn't.

nT = min(diff(cropT,1,2));

dir_now = pwd;

cd(dataPath)
if exist('o3','dir')
    error('ozone!')
end

system('mkdir o3 && mv *.dat o3')

cd o3

dl = dir;
dl = dl(3:end);

for idx = 1:length(dl)
    fn = dl(idx).name(1:end-4);
    for jdx = 1:size(cropT,1)
        com = sprintf('dd if=%s.dat of=../%s%i.dat ibs=8 count=%i skip=%i', fn, fn, jdx, nT, cropT(jdx,1)-1);
        [~, result] = system(com);
    end
end

cd ../
ls -l
cd(dir_now)

pwd