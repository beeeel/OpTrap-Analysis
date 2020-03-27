%% Show unwrapped result with errorbars
% Show errors from standard deviation of stationary deformation. 

SubPlots = true;
SvPng = true;
SvFig = true;

CellType = 'LS174T';
Set = 'normoxia';
Fh = figure(9);
Nums = [2 9 11 14 16 19];
for Num = 1:length(Nums)
    if Subplots
        subplot(length(Nums),1,Num)
        cla
    end
    [~, ~, ~] = PlotUnwrapErrors(CellType, Set, num2str(Nums(Num)),false,false);
end

if Subplots
    %% Fixup if subplots
    Axs = ones(length(Fh.Children));
    for Ax = Fh.Children'
        This = find(Fh.Children == Ax);
        if This ~= 6
            Ax.Title.String = Ax.Title.String(end);
        end
        if This ~=3
            Ax.YLabel.String = '';
        end
        if This ~= 1
            Ax.XLabel.FontSize = 8;
        end
    end
end

%% Save
if SvPng
    saveas(gcf,['~/png/EQ/D_error_std_multiplot_' strjoin({CellType, Set},'_') '.png'])
end
if SvFig
    saveas(gcf,['~/fig/EQ/D_error_std_multiplot_' strjoin({CellType, Set},'_') '.fig'])
end