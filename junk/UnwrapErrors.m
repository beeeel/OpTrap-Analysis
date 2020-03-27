%% Show unwrapped result with errorbars
% An attempt at characterising errors from unwrap analysis compared to the
% "gold standard" errors from fitting to simulated data



CellType = 'LS174T';
Set = 'normoxia';
Fh = figure(9);
Nums = [2 9 11 14 16 19];
for Num = 1:length(Nums)
    subplot(length(Nums),1,Num)
    cla
    [~, ~, ~] = PlotUnwrapErrors(CellType, Set, num2str(Nums(Num)),false,false);
end


