function SaveFigPng(FileName, Dir, SaveFig, SavePng)
%% SaveFigPng Save a figure as .fig and .png into hardcoded directories
% SaveFigPng(FileName, Dir, SaveFig, SavePng) 
% FileName - char. Dir - char. SaveFig - logical. SavePng - logical.

if ~strcmp(Dir(end),'/')
    Dir = [Dir '/'];
end

FigDir = ['~/fig/' Dir];
PngDir = ['~/png/' Dir];

if SaveFig
    saveas(gcf, [FigDir FileName '.fig'])
end

if SavePng
    saveas(gcf, [PngDir FileName '.png'])
end
end