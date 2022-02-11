function eventAnimator(cdEvents, ROI, t0, step, window, nSteps, outputPath)

if ~exist('ROI','var')
    ROI = [0 0 640 480];
end
if ~exist('t0','var')
    t0 = 0;
end
if ~exist('step','var')
    step = 0.1;
end
if ~exist('window','var')
    window = 0.1;
end
if ~exist('nSteps','var')
    nSteps = 100;
end
if ~exist('outputPath','var')
    outputPath = '';
elseif exist(outputPath,'file')
    error('Not overwriting %s',outputPath)
else
    vid = VideoWriter(outputPath,'Motion JPEG AVI');
    open(vid);
end

% All event indices to animate
idxs = [find(cdEvents.ts*1e-6>t0,1,'first') find(cdEvents.ts*1e-6<t0 + nSteps*step+window,1,'last')]; 
figure(4)
clf
colormap gray
for idx = 1:nSteps
    im = zeros(480, 640);
    % Calculate the indices 
    jdxs = ( (cdEvents.ts*1e-6 > t0 + (idx-1)*step) & (cdEvents.ts*1e-6 < t0 + (idx-1)*step+window));
    i1 = find(jdxs, 1, 'first');
    i2 = find(jdxs, 1, 'last');
    
    % Make the image
    for jdx = i1:i2
        y = 1+cdEvents.y(jdx); 
        x = 1+cdEvents.x(jdx);
        im(y, x) = im(y,x) + cdEvents.p(jdx);
    end
    % Find the centroid of the events
    if ~isempty([i1 i2])
        if i1 ~= i2
            [t, c] = bufferEvents_2(...
                cdEvents.ts(jdxs), ...
                cdEvents.x(jdxs), ...
                cdEvents.y(jdxs), ...
                sum(jdxs), ROI);
            % This index is always correct if there are events are in
            % the buffer.
            if ~isempty(t)
                t = t(round(length(t)/2));
                c = c(:,round(length(t)/2));
            end
        else
            t = [];
            c = [];
        end
    else
        t = [];
        c = [];
    end

%     if idx == 10
%         error('aaa')
%     end
    subplot(20, 10, [1 10])
    cla
    hold on
    title(sprintf('Events #%i to #%i',i1, i2))
    ylabel('#')
    xlabel('Time (s)')
    plot(cdEvents.ts(idxs)./1e6, idxs,'k', 'LineWidth', 2)
    plot(cdEvents.ts([i1 i2])./1e6, [i1 i2], 'y', 'LineWidth', 2)
    if ~isempty(t)
        plot(t./1e6, (i1+i2)/2, 'rx','MarkerSize',6,'LineWidth',3)
    end
    subplot(20, 10, [21 200])
    imagesc(im, [-1 1]*3);
    hold on
    plot(ROI(1).*ones(1,5)+[0 1 1 0 0].*ROI(3), ROI(2).*ones(1,5)+[0 0 1 1 0].*ROI(4), 'r-', 'LineWidth',2)
    if ~isempty(c)
        plot(c(1), c(2), '.r', 'LineWidth', 3, 'MarkerSize', 9)%, 'Color', [1 1 1 1] * 0.3)
    end
    
    colorbar
    axis equal
    if ~any(ROI([3,4]) == [640 480])
        xlim(ROI(1) + [-0.5 1.5]*ROI(3))
        ylim(ROI(2) + [-0.5 1.5]*ROI(4))
    end
    title(sprintf('%.3fs to %.3fs, window length %.3gs',1e-6*cdEvents.ts(i1),1e-6*cdEvents.ts(i2), window))
%     fprintf('%.1fs to %.1fs, %g events\n',1e-6*cdEvents.ts((idx-1)*step + 1),1e-6*cdEvents.ts((idx-1)*step + step+2*overlap), step+2*overlap)
    if ~isempty(outputPath)
        cf = getframe(gcf);
        writeVideo(vid, cf);
    else
        pause(0.015)
    end
end

if ~isempty(outputPath)
    close(vid);
end
   