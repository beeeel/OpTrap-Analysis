clear all; close all; clc;

saveVar = 1;

%-----------user prompt to choose folder to read files from------------------
[chosenFiles, chosenpath] = uigetfile('*.avi', 'Select a video','MultiSelect','on');
fileNames = fullfile(chosenpath, chosenFiles);


%-----------user prompt to choose folder to save files to------------------
[savePath] = uigetdir('Select a folder to save the output data in');

for idx = 1:length(fileNames(1,:));
    
    if length(fileNames(1,:))>1
        fileName = char(fileNames(idx));
    elseif length(fileNames(1,:))==1
        fileName = fileNames;
    else
        return
    end
    
    %PostProcess Video
    CellDeformInfo = PostProcessCellDeform(fileName);
    
    %-----------save data to file in folder specified--------
    if saveVar == 1
        [~,fileSaveName,~] = fileparts(fileName);
        
        fileSaveName = strcat(fileSaveName,'Processed');
        fullSaveName = fullfile(savePath,fileSaveName);
        save(fullSaveName,'CellDeformInfo');
    end
    
end