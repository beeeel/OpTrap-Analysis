dataDir = 'E:/Will/data/2021_06_30/';
addpath E:/Will/OpTrap-Analysis/helpers/
addpath E:/Will/OpTrap-Analysis/beads/
close all

lastLength = 2;
while true
    dirList = dir(dataDir);
    if length(dirList) ~= lastLength
        nFiles = length(dirList) - lastLength;
        if nFiles < 0
            error('You deleted something and were too lazy to code this properly')
        end
        
        for fIdx = nFiles%:-1:1
            idx = length(dirList)-fIdx+1;
            if dirList(idx).isdir
                close
                fprintf('About to load %s \n', dirList(idx).name)
                data = struct();
                data.dirPath = ([dataDir dirList(idx).name]);
                data.fName = dirList(idx).name;
                data.mPerPx = 7e-8;
                data.opts.cropT = [];
                data.opts.pOrder = 1;
                data = bead_loadData(data, false);
                %data = bead_preProcessCentres(data);
                
                bead_plotRawData(data);
                %bead_normMSD_polyfit( data, 'a', 1);
            end
        end
        
        
        lastLength = length(dirList);
    end
    
    pause(5)
    
end


