function bead_plotCompareCentres(data, varargin)
%% plotCompareCentres(data, [step])
% Shows images from Imstack, centres from [xCentres yCentres]
% overlaid along with calculated centres using normalize-square CofM
% algorithm in MATLAB. Step is used to select correct frames from centres
% arrays.

Imstack = data.Imstack;
xCentres = data.raw.xCentresPx;
yCentres = data.raw.yCentresPx;
timeVecMs = data.raw.timeVecMs;

if nargin == 1
    step = length(xCentres)./size(Imstack{1},1);
elseif nargin == 2
    step = varargin{1};
else
    error('Nargins should be 1 or 2!')
end
% Get two centres arrays
imCentres = imCentreOfMass(cat(3,Imstack{1}{:,1}));
imCentres2 = imCentreOfMass(cat(3,Imstack{1}{:,1}),'dark-norm-square');
imCentres = [imCentres; imCentres2];
xyCentres = cat(3,xCentres(:,1:step:end), yCentres(:,1:step:end));
times = timeVecMs(1:step:end);
% Create a UIFigure with axes holding plots and plot the first
% image
fh = uifigure('Name','Image with both calculated centres',...
    'Position',[680 160 980 800]);
colormap(fh, 'gray');
ax = axes(fh);
hold(ax,'on')
changeIm(struct('Value',1), Imstack, ax, imCentres, xyCentres, times, data.raw.suffixes);
% Create a slider for user control of which image is showing
sld = uislider(fh, 'Position', [100, 50, 600, 40], ...
    'ValueChangedFcn', @(sld, event) changeIm(sld, Imstack, ax, imCentres, xyCentres, times, data.raw.suffixes),...
    'Limits', [1 size(Imstack{1},1)], 'MinorTicks', 1:size(Imstack{1},1));
input('Enter to close figure window')
close(fh)
end

function changeIm(sld, ims, ax, imCentres, xyCentres, times, suffixes) 
fr = round(sld.Value);
cla(ax);
imagesc(ax, ims{1}{fr,1});
plot(ax, imCentres(1, fr), imCentres(2, fr), 'yx','MarkerSize',6);
plot(ax, imCentres(3, fr), imCentres(4, fr), 'gx','MarkerSize',6);
legCell = {'MATLAB simple centre','MATLAB dark centre'};
for idx = 1:size(xyCentres,1)
    plot(ax, xyCentres(idx, fr, 1), xyCentres(idx, fr, 2),'.','MarkerSize',12);
    legCell{end+1} = ['Live: ' suffixes{idx}];
end
axis(ax, 'image');
legend(ax, legCell)
diffCentres = imCentres - xyCentres(1, 1:length(imCentres),:);
title(ax,{['Frame ' num2str(fr) ' x_{diff} = ' num2str(diffCentres(1,fr)) ...
    ' y_{diff} = ' num2str(diffCentres(2,fr))], ['t = ' num2str(times(fr)/1e3) 's']});
end