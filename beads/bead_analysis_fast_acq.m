%% Load real-time processed bead data from B09
filePaths = {'~/Documents/data/OpTrap/bead_analysis/2020_11_16/500nm_50pc_laser_1/',...
    '~/Documents/data/OpTrap/bead_analysis/2020_11_16/500nm_50pc_laser_2/',...
    '~/Documents/data/OpTrap/bead_analysis/2020_11_12_2um/1/',...
    '~/Documents/data/OpTrap/bead_analysis/2020_11_12_2um/2/',...
    '~/Documents/data/OpTrap/bead_analysis/2020_11_12_2um/3/',...
    '~/Documents/data/OpTrap/bead_analysis/2020_11_12_2um/4'};
mPerPx = 0.07e-6;
showStack = false;
cropTs = {[1 288776], [8e4 228e3],[1 1e6],[1 139e3],[1 1e6],[1 1e6]};
for fileIdx = 4%1:length(filePaths)
    xCentres = byteStreamToDouble([filePaths{fileIdx} 'XCentres.dat']);
    yCentres = byteStreamToDouble([filePaths{fileIdx} 'YCentres.dat']);
    timeVec  = byteStreamToDouble([filePaths{fileIdx} 'Times.dat']);
    Imstack  = bfopen([filePaths{fileIdx} 'images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    metadata = fileread([filePaths{fileIdx} 'images_and_metadata/images_and_metadata_MMStack_Default_metadata.txt']);
    metadata = jsondecode(metadata);
    fName = strsplit(filePaths{fileIdx},'/');
    fName = fName{end-1};
    %%
    cropT = cropTs{fileIdx};
    xCentresM = xCentres(cropT(1):cropT(2)) .* mPerPx;
    yCentresM = yCentres(cropT(1):cropT(2)) .* mPerPx;
    xStiff = calcStiffness(xCentresM);
    yStiff = calcStiffness(yCentresM);
    
    imCentres = imCentreOfMass(cat(3,Imstack{1}{:,1}));
    
    fh = figure(fileIdx);
    fh.Name = fName;
    clf
    subplot(3,1,1)
    hold on
    histogram(xCentresM.*1e6)
    histogram(yCentresM.*1e6)
    xlabel('Centre position (\mu m)')
    ylabel('Bin count')
    title(['Histogram of centres, trap stiffness kx = ' num2str(xStiff./1e-6) ' pN/\mu m, ky = ' num2str(yStiff./1e-6) ' pN/\mu m'])
    legend('X','Y')
    subplot(3,1,2)
    plot(xCentresM.*1e6,yCentresM.*1e6,'.')
    title(['Scatterplot of centres, ' num2str(diff(cropT)+1) ' frames'])
    xlabel('X Centre position (\mu m)')
    ylabel('Y Centre position (\mu m)')
    axis equal
    subplot(3,2,5)
    plot(1e-3.*timeVec(cropT(1):cropT(2)), xCentresM,'.')
    xlabel('Time (s)')
    ylabel('X centre (\mu m)')
    title('X centre position time trace')
    subplot(3,2,6)
    plot(1e-3.*timeVec(cropT(1):cropT(2)), yCentresM,'.')
    xlabel('Time (s)')
    ylabel('Y centre (\mu m)')
    title('Y centre position time trace')
    drawnow
    if showStack
        disp('MATLAB is locked until you close imagej')
        system(['imagej ' filePaths{fileIdx} 'images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    end
end

%% Show calculated centres with image
frames = 1:15;

imCentres = imCentreOfMass(cat(3,Imstack{1}{frames,1}));
figure(4)
fh.Name = fName;
clf
hold on
plot(0.07.*xCentres(1:1e4:max(frames)*1e4),0.07.*yCentres(1:1e4:max(frames)*1e4),'k.')
plot(0.07.*imCentres(1,:), 0.07.*imCentres(2,:), 'rx')
axis equal

figure(5)
fh.Name = fName;
for frame = 1%:16
clf
hold on
imagesc(Imstack{1}{frame,1})
plot(imCentres(1,frame), imCentres(2,frame),'x')
plot(xCentres((frame-1).*1e4+1)+0.5,yCentres((frame-1).*1e4+1)+0.5,'.')
legend('MATLAB measured centre','ImageJ measured centre')
title(num2str(frame))
pause(0.5)
end
%% Methods development - thresholding
frame = 1;
xLine = 34;
im = Imstack{1}{frame,1};
im(im<2*mean(im,'all')) = 0;
figure(1)
imshowpair(Imstack{1}{frame,1},im)
figure(2)
plot(1:size(Imstack{1}{1,1},2),im(xLine,:))