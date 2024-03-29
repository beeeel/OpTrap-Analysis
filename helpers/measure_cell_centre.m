function cCentre = measure_cell_centre(im, dirPath, roi)
% Measure the centre position of a cell in an image, and store the result
% in a text file in the given directory. cCentre contains [Xc, Yc].

if ~strcmp(dirPath(end), '/')
    dirPath(end+1) = '/';
end

cFile = [dirPath '/cell_centre.txt'];
if exist(cFile, 'file')
    warning('Cell centre measurement already exists!')
    i = input('Overwrite previous measurement? (y/[N])', 's');
    if ~strcmp(i, 'y')
        warning('Received input %s, aborting', i)
        return
    end
end
fid = fopen(cFile, 'w');


%%
fh = figure('Units','normalized','OuterPosition',[0 0 1 1]);
clf
imagesc(im, prctile(reshape(im, [], 1), [5 95]))
hold on
set(gca,'FontSize',16)
axis image
axis off
plot(roi(1) * ones(1,5) + roi(3) * [0 1 1 0 0], roi(2) * ones(1,5) + roi(4) * [0 0 1 1 0], 'r-','LineWidth',2)

title({'Drawn an ellipse on the cell,' 'matching the curvature of the cell where the bead touches'});
hUI = uicontrol(fh, 'Style', 'togglebutton', 'String', 'Done','Position',[20 20 180 120], ...
    'BackgroundColor',[1 1 1]/4, 'FontSize', 36, 'ForegroundColor',[1 1 1], ...
    'FontWeight', 'bold');
ellipse = drawellipse('Color','r');
title('Press ''Done'' when you are finished')

while ~hUI.Value
    pause(1)
end

cCentre = ellipse.Center;
title('Measured!')

close(fh)

fprintf(fid, '%g  %g', cCentre);
fclose(fid);
end