%% Manually measure the size of a cell

% Load file
dirPath = '~/Documents/data/OpTrap/2020_11_19/hela_010mms-1_60umdeep/';
fileName = 'hela_010mms-1_60umdeep_MMStack_Default.ome.tif';
Imstack = bfopen([dirPath fileName]);
[imH, imW] = size(Imstack{1}{1,1});
metadata = fileread([dirPath fileName(1:end-8) '_metadata.txt']);
metadata = jsondecode(metadata);
%% Manually choose two frames from the beginning of the drag and the end
frames = [1, 680];
frameTimes = zeros(size(frames));
for idx = 1:length(frames)
    frameTimes(idx) = metadata.(['FrameKey_' num2str(frames(idx)-1) '_0_0']).ElapsedTime_ms;
end
%%
figure(1)
clf
imagesc([1 imW].*0.07, [1 imH].*0.07, Imstack{1}{frames(1),1})
xlabel('\mu m')
ylabel('\mu m')
axis image
roi = drawpoint;
close
%%
lines(:,1) = Imstack{1}{frames(1),1}(:,round(roi.Position(1)./0.07));
lines(:,2) = Imstack{1}{frames(2),1}(:,round(roi.Position(1)./0.07));
figure(2)
clf
for plt = 1:2
    subplot(2,2,plt)
    hold on
    imagesc([1 imW].*0.07, [1 imH].*0.07, Imstack{1}{frames(plt),1})
    axis image
    xlabel('\mu m')
    ylabel('\mu m')
    plot(roi.Position(1).*[1 1], [1 imH].*0.07,'r')
    title(['Frame ' num2str(frames(plt)) ' T_{drag} = ' num2str(frameTimes(plt)./1e3 - 1) 's'])
end
subplot(2,1,2)
hold on
plot((1:imH).*0.07, lines(:,1),'k')
plot((1:imH).*0.07, lines(:,2),'r')
xlabel('Position (\mu m)')
ylabel('Pixel value')
legend(['Frame ' num2str(frames(1))],['Frame ' num2str(frames(2))])