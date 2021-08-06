function cCentre = measure_cell_centre(im, dirPath)
% Measure the centre position of a cell in an image, and store the result
% in a text file in the given directory

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
    fid = fopen(cFile, 'w');
end

%%
fh = figure('Units','normalized','OuterPosition',[0 0 1 1]);
clf
imagesc(im, prctile(reshape(im, [], 1), [5 95]))
set(gca,'FontSize',16)
axis image
axis off
title({'Drawn an ellipse on the cell,' 'matching the curvature of the cell where the bead touches'});
hUI = uicontrol(fh, 'Style', 'togglebutton', 'String', 'Done','Position',[20 20 180 120], ...
    'BackgroundColor',[1 1 1]/4, 'FontSize', 36, 'ForegroundColor',[1 1 1], ...
    'FontWeight', 'bold');
ellipse = drawellipse;
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