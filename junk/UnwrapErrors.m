%% Show unwrapped result with errorbars
% Show errors from standard deviation of stationary deformation. 

SubPlots = true;
SvPng = true;
SvFig = true;

CellType = 'HL60';
Set = 'with_drugs';
Fh = figure(9);
Nums = [4,5,6,8,9,11,13,18,20];
for Num = 1:length(Nums)
    if SubPlots
        subplot(ceil(length(Nums)/2),2,Num)
    end
    cla
    [~, ~, ~] = PlotUnwrapErrors(CellType, Set, num2str(Nums(Num)),false,false);
    SaveFigAndPng(CellType, Set, SvPng, SvFig, Num, 0)
end

if SubPlots
    %% Fixup if subplot
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
    SaveFigAndPng(CellType, Set, SvPng,SvFig, 0, false)
end

%%
function SaveFigAndPng(CellType, Set, SvPng, SvFig, Num, Scaled)
%% Save
if Scaled
    SC = '_sc';
else
    SC = '';
end

if ~Num
    Fname = ['_multiplot' SC '_std_D'];
else
    Fname = ['_' num2str(Num) SC '_std_D'];
end

if SvPng
    saveas(gcf,['~/png/EQ/D_' strjoin({CellType, Set},'_') Fname '.png'])
end
if SvFig
    saveas(gcf,['~/fig/EQ/D_' strjoin({CellType, Set},'_') Fname '.fig'])
end
end