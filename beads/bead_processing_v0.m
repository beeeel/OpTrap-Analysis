dirList = dir;
dirList = dirList([dirList.isdir]);
dirList = dirList(3:end);
mPerPx = 0.07e-6;
showStack = false;
doPlots = true;
cropTs = {[1 5e5], [1 5e5], [1 5e5], [1e5 3.2e5], [1.4e5 4.4e5], [1.7e5 5e5], [1 0.9e5], [1.6e5 2.8e5] };
laserPowers = 30:5:60;
stiffXY = zeros(2,length(laserPowers));
checkCropTs(cropTs, dirList);
for fileIdx = 1:length(dirList)
    dirPath = [dirList(fileIdx).folder '/' dirList(fileIdx).name];
    xCentres = byteStreamToDouble([dirPath '/XCentres.dat']);
    yCentres = byteStreamToDouble([dirPath '/YCentres.dat']);
    timeVec  = byteStreamToDouble([dirPath '/Times.dat']);
    Imstack  = bfopen([dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    metadata = fileread([dirPath '/images_and_metadata/images_and_metadata_MMStack_Default_metadata.txt']);
    metadata = jsondecode(metadata);
    fName = dirList(fileIdx).name;
    %%
    cropT = cropTs{fileIdx};
    xCentresM = xCentres(cropT(1):cropT(2)) .* mPerPx;
    yCentresM = yCentres(cropT(1):cropT(2)) .* mPerPx;
    xStiff = calcStiffness(xCentresM);
    yStiff = calcStiffness(yCentresM);
    stiffXY(:, fileIdx) = [xStiff, yStiff];
    
    imCentres = imCentreOfMass(cat(3,Imstack{1}{:,1}));
    
    if doPlots
        fh = figure(fileIdx); %#ok<*UNRCH>
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
    end
    if showStack
        disp('MATLAB is locked until you close imagej')
        system(['imagej ' dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    end
end
%% Plot stiffnesses of all sets
figure(fileIdx+1)
clf
subplot(2,1,1)
hold on
plot(laserPowers,1e6*stiffXY(1,1:length(laserPowers)))
plot(laserPowers,1e6*stiffXY(2,1:length(laserPowers)))
title('Stiffness varying with laser power setting')
legend('X stiffness', 'Y stiffness')
xlabel('Laser power setting (%)')
ylabel('Trap stiffness (pN/\mu m)')
subplot(2,1,2)
bar(laserPowers, stiffXY(1,1:length(laserPowers))./stiffXY(2,1:length(laserPowers)))
title('Ratio $\frac{k_x}{k_y}$','Interpreter','latex','Fontsize',20)
xlabel('Laser power setting (%)')
ylabel('Ratio')
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

function checkCropTs(cell, struct)
if length(cell) ~= length(struct)
    str = 'cropTs = {';
    str = [str repmat('[1 5e5], ', 1, length(struct))];
    str = [str(1:end-2) ' };'];
    disp(str);
    error('CropTs is empty, copy the above line into the script');
end
end