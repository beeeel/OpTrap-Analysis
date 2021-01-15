%% Find/segment the cell assisted by background subtraction
CellType = 'HL60';
Set = 'normoxia';
Num = '2';

global Imstack info meta
LoadImstackInfoMeta(CellType,Set,Num,true)
%%
ImsBG = double(mean(cat(3,Imstack{1}{1:10,1}),3));

ImsFilt = double(imgaussfilt(cat(3,Imstack{1}{:,1}),1));
ImsDiff = (ImsFilt - ImsBG).^2;
%%
figure(1)
clf
subplot(2,1,1)
imagesc(Imstack{1}{1,1})
axis image off
subplot(2,1,2)
imagesc(ImsDiff(:,:,200))
axis image off