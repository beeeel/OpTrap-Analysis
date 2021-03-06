%% Show unwrapped result with errorbars
% Show errors from standard deviation of stationary deformation. 
SubPlots = true;
SvMulti = true;
SvPng = false;
SvFig = false;

CellType = 'HL60';
Set = 'normoxia';
Nums = [1:20];
% Nums = [4,5,6,8,9,11,13,18,20]; % HL60 drugs nums
%%
Fh = figure(9);
clf
for Num = 1:length(Nums)
    if SubPlots
        subplot(ceil(length(Nums)/4),4,Num)
    end
    cla
    PlotUnwrapErrors(CellType, Set, num2str(Nums(Num)),false,false);
    if ~SubPlots
        SaveFigAndPng(CellType, Set, SvPng, SvFig, Num, 0)
    end
end

TtlNum = 19;
% XLNum = 12;
YLabNum = 12;
if SubPlots
    %% Fixup if subplot
    Axs = ones(length(Fh.Children));
    for Ax = Fh.Children'
        ThisAxNum = find(Fh.Children == Ax);
        if ThisAxNum ~= TtlNum
            Ax.Title.String = Ax.Title.String(end);
        end
        if ThisAxNum ~=YLabNum
            Ax.YLabel.String = '';
        end
        Ax.XLabel.FontSize = 14;
%         if This ~= 1
%             Ax.XLabel.FontSize = 8;
%         end
    end
    SaveFigAndPng(CellType, Set, SvMulti,SvMulti, 0, false)
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
