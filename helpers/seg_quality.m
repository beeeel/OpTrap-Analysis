%% Notes:
% This needs to be turned into a few functions - a loader function, a
% figure function with proper functionality (overlays, 2-6 plots, choice on
% what the plots are)
%% Load data
CellType = 'LS174T';
Set = 'normoxia';
Num = '11';
datadir = '/home/will/Documents/data/OpTrap/';
infosdir = [datadir 'infos/'];
if isempty(whos('Imstack')); Imstack = {{0,''}}; end % Create an empty

if strcmp(CellType,'hela')
    load([infosdir 'info-hela_ctrl_s_020_tr_70_' Num '.mat'])
    imfile = [datadir '0610/Deformation/hela_ctrl_s_020_tr_70_' Num '_MMStack.ome.tif']
    if ~strcmp(imfile, Imstack{1}{1,2}(1:length(imfile)))
        Imstack = load_imstack('0610/Deformation','ctrl','020','70',0,Num);
    end
elseif strcmp(Set,'Jenna')
    error('you need to update this branch of the load operation')
    load(strcat('/home/ppxwh2/Documents/data/OpTrap/infos/info_Jenna_test_no_erode',Num,'.mat'))
    imfile = ['/home/ppxwh2/Documents/data/OpTrap/1119_Jenna/191119_thp1_ctrl_s_010_tr_50_',Num,'_MMStack.ome.tif']
    if ~strcmp(imfile, Imstack{1}{1,2}(1:length(imfile)))
        Imstack = bfopen(imfile);
    end
elseif strcmp(CellType,'HL60')
    if strcmp(Set,'with_drugs')
        matName = ['_190717_HL60_' Num '_0.020mm-1_1'];
    elseif strcmp(Set,'normoxia')
        matName = ['_100717_' Num '_' CellType '_1']; 
    end
    load(['/home/ppxwh2/Documents/data/OpTrap/infos/info_' CellType '_' Set matName '.mat'])
    imfile = ['/home/ppxwh2/Documents/data/OpTrap/2017_10_movies-from-aishah/'...
        CellType '/' CellType '_' Set '/' matName(2:end) '.avi']
    if ~strcmp(imfile, Imstack{1}{1,2}(1:length(imfile)))
        Imstack = avi_to_imstack(imfile)
    end
elseif strcmp(CellType,'LS174T')
    imfile = ['200717_' Num '_LS174T_' Set '_1.avi'];
    if strcmp(Set,'hypoxia'); imfile(2) = '1'; end
    load([infosdir 'info_seg_LS174T_' Set '_' imfile(1:end-4) '.mat'])
    Imstack = avi_to_imstack([datadir '2017_10_movies-from-aishah/LS174T/' imfile]);
end

%% Segmentation video
% This can also save a gif or an avi movie
% Opens a figure window with an overlay of the segmentation atop the image,
% and either 4 or 6 plots (see n_plots). The 4 plots show Taylor Parameter,
% Area and Centroids from regionprops, and either orientation from
% regionprops or radius from find_cell. The extra 2 are either the crop box
% co-ordinates from find_cell and the success metrics from both modules, or
% the results of the ellipse fitting (Taylor Parameter and orientation).
%
% The mask overlay can either be the initial circular mask from
% find_circle, or the final mask from segment_cell. Atop this can overlay
% the fitted ellipse from regionprops or ellipseDetection.
%
% The time to pause on each frame can be set, as can save options to create
% a .gif or .avi from the slideshow. If flythrough is set to 0, you need to
% press enter on the command window to advance each frame.

show_fit = 'unwrap';   % 'regionProps' or 'ellipseDetection' or 'unwrap' or 'none' or 'linemax' - source of overlay on top of mask (linemax is centre only)
show_mask = 'none';      % 'initial' or 'segment' or 'none' - Segmented mask from seg_cell, or initial circular mask from find_cell
n_plots = 0;                % Number of plots - 6 includes ellipse fitting results
pt_mode = 'data';           % Analysis or data or unwrap - do you want to look at the data, or analyse why it isn't working, or just show unwrapped data
frs =50;                 % Frames to display

p_time = 0.25;              % Time to pause on each frame when showing as movie
makevid = 0;                % Set to 1 to make animated gif or 2 to make avi
flythrough = 1;             % Play as movie instead of requiring user input to step through
run_unwrap = 0;             % Run unwrap cell before displaying data (refreshes content of unwrapped, fits, Ia)

% Font sizes for axes and titles
XFontSize = 12;
YFontSize = 12;
TFontSize = 12;

if isempty(whos('offset')); run_unwrap = 1; end
if run_unwrap && meta.find_cell_v; [u_fits, unwrapped, Ia, FitEqn, offset] = unwrap_cell_v2(Imstack, [info.centres] , [info.radius],'sc_up',1.8,'ifNaN','mean','sc_down',0.35); 
elseif run_unwrap; [u_fits, unwrapped, Ia, FitEqn, offset, FitErrs] = unwrap_cell_v2(Imstack, [info.mCentres] , repmat(100,1,size(Imstack{1},1)),'sc_up',1.8,'ifNaN','mean','sc_down',0.35); end %#ok<UNRCH>

if n_plots == 4; sbplt = [4 2 4]; 
elseif n_plots == 6; sbplt = [6 2 6]; end

filename = strsplit(info(1).filepath,{'/','.','_'});
fname = strjoin({'seg',filename{9:15}},'_'); % Filename base for saving - seg_(original filename)
h = figure(13);


centroids = [info.Centroid];
if run_unwrap
    for fr = 1:length(Imstack{1})
        info(fr).uMajorAxisLength = u_fits(1,fr);
        info (fr).uMinorAxisLength = u_fits(2,fr);
        info(fr).uOrientation = u_fits(3,fr);
        info(fr).uTaylorParameter = (u_fits(1,fr) - u_fits(2,fr))./(u_fits(1,fr) + u_fits(2,fr));
    end
end

if meta.seg_cell_v == 5; fits = [info.ellipse_fits]; end
theta = 0:0.01:2*pi;
if strcmp(pt_mode,'analysis'); crops = [info.crop]; end
if strcmp(show_mask,'initial')
    [rr, cc] = meshgrid(1:size(Imstack{1}{1,1},2),...
        1:size(Imstack{1}{1,1},1));
    centres = [info.centres];
elseif strcmp(show_fit,'unwrap') && meta.line_maxima_v
    centres = [info.mCentres];
end
Xdata = 1:size(info,2);

if makevid == 2; v = VideoWriter(strcat(fname,'.avi')); v.FrameRate = 4; open(v); end
for frame = frs
    if length(whos('unwrapped','Ia','fiteqn')) == 3
        subplot(2,2,2)
        imagesc(1:360,(1:size(unwrapped,1))*0.07,unwrapped(:,:,frame))
        hold on
        plot(0.07 * Ia(:,:,frame),'r.')
        plot(0.07 * FitEqn(info(frame).uMajorAxisLength, info(frame).uMinorAxisLength, ...
            info(frame).uOrientation, 1:360),'k:','LineWidth',3)
        hold off
        title('Unwrapped cell with fitting data and Centroid fitted data','FontSize',TFontSize)
        xlabel('CW angle from +x (degrees)','FontSize',XFontSize)
        ylabel('Radius (\mum)','FontSize',YFontSize)
        legend('Column max','Fitted ellipse')
        % add grid, change xticks
        ax = gca;
        ax.XTick = 0:90:360;
        grid on
        subplot(2,2,1)
    else
        subplot(2,1,1)
    end
    hold off
    if strcmp(show_mask,'segment')
        imshowpair(Imstack{1}{frame,1},info(frame).mask)
        hold on
    elseif strcmp(show_mask,'initial')
        % Note: This isn't the same mask as used by segment_cell: the value
        % for sc_up is not saved, so this just uses 1.25 (default)
        imshowpair(Imstack{1}{frame,1},((cc-centres(2,frame)).^2 + ...
            (rr-centres(1,frame)).^2 <= (info(frame).radius*1.25).^2));
        hold on
        plot(centres(1,frame),centres(2,frame),'y.','MarkerSize',16)
    elseif strcmp(show_mask,'none')
        aa = imagesc(Imstack{1}{frame,1});
        colormap(aa.Parent,'gray')
        hold on
    end
    axis image off
    if size(Imstack{1}{1,1},1) == 1080
        ylim([270, 810]), xlim([480, 1440])
    end
    title(['Segmentation and fitting: frame ' num2str(frame)],'FontSize',TFontSize+4)
    % Draw the ellipse on - for an ellipse with centre (x0, y0), semi-axis
    % lengths (a,b), oriented at an angle phi above horizontal, the
    % cartesian equations from the polar are as follows:
    % x = a cos(theta) cos(phi) - b sin(theta) sin(phi) + x0
    % y = a cos(theta) sin(phi) + b sin(theta) cos(phi) + y0
    % Which is translated into indexed variables below
    
    if strcmp(show_fit, 'regionProps')
        % Using regionprops fitting
        plot(0.5 * info(frame).MajorAxisLength .* cos(theta) .* cos(info(frame).Orientation) ...
            - 0.5 * info(frame).MinorAxisLength .* sin(theta) .* sin(info(frame).Orientation)...
            + info(frame).Centroid(1),... % x values end here
            0.5 * info(frame).MajorAxisLength .* cos(theta) .* sin(info(frame).Orientation) ...
            + 0.5 * info(frame).MinorAxisLength .* sin(theta) .* cos(info(frame).Orientation)...
            + info(frame).Centroid(2),'k--','LineWidth',2)
        plot(centroids(frame*2-1),centroids(frame*2),'kx')
    elseif strcmp(show_fit, 'ellipseDetection')
        % Using ellipseDetection fitting
        plot( info(frame).fits(3) .* cos(theta) .* cos(info(frame).fits(5)) ...
            - info(frame).fits(4) .* sin(theta) .* sin(info(frame).fits(5))...
            + info(frame).fits(1),... % x values end here
            info(frame).fits(3) .* cos(theta) .* sin(info(frame).fits(5)) ...
            + info(frame).fits(4) .* sin(theta) .* cos(info(frame).fits(5))...
            + info(frame).fits(2),'y','LineWidth',2)
        plot(fits(1,frame),fits(2,frame),'kx')
    elseif strcmp(show_fit, 'unwrap')
        % Using unwrap cell fitting
        plot(info(frame).uMajorAxisLength .* cos(theta) .* cos(info(frame).uOrientation) ...
            - info(frame).uMinorAxisLength .* sin(theta) .* sin(-info(frame).uOrientation)...
            + centres(1,frame) + offset(2,frame),... % x values end here
            info(frame).uMajorAxisLength .* cos(theta) .* sin(-info(frame).uOrientation) ...
            + info(frame).uMinorAxisLength .* sin(theta) .* cos(info(frame).uOrientation)...
            + centres(2,frame) + offset(3,frame),'k--','LineWidth',2)
        plot(centres(1,frame) + offset(2,frame),centres(2,frame) + offset(3,frame),'kx')
    elseif strcmp(show_fit,'linemax')
        plot(Centres(1,frame),Centres(2,frame),'kx')
    end
    if n_plots
        subplot(sbplt(1), sbplt(2), sbplt(3) + 1)
        trackPlot(Xdata,[info.TaylorParameter],frame)
        title('Taylor Parameter [regionprops]','FontSize',TFontSize)
        YMax = 1.1*max([info.TaylorParameter]);
        if isnan(YMax); YMax = 1; end
        ylim([0, YMax])
        
        subplot(sbplt(1), sbplt(2), sbplt(3) + 2)
        trackPlot(Xdata,(0.07^2)*[info.Area],frame)
        title('Area [regionprops]','FontSize',TFontSize), ylabel('\mu m^2','FontSize',YFontSize)
        YMax = 1.1*max([info.Area])*(0.07^2);
        if isnan(YMax); YMax = 1; end
        ylim([0, YMax])
        
        subplot(sbplt(1), sbplt(2), sbplt(3) + 3)
        hold off
        plot(centroids(1:2:end)*0.07)
        hold on
        plot(centroids(2:2:end)*0.07)
        plot(frame,centroids(2*frame - 1)*0.07,'rx')
        plot(frame,centroids(2*frame)*0.07, 'rx')
        ylim([0 length(Imstack{1}{1,1})]*0.07), ylabel('\mu m','FontSize',YFontSize)
        title('Centroids of mask [regionprops]','FontSize',TFontSize)
        if n_plots == 4; xlabel('frame number'); end
        
        subplot(sbplt(1), sbplt(2), sbplt(3) + 4)
        if strcmp(pt_mode,'data')
            trackPlot(Xdata,[info.Orientation],frame,'x')
            title('Angle of long axis above horizontal [regionprops]','FontSize',TFontSize)
            ylabel('Angle (^{o})','FontSize',YFontSize)
            ylim([-90 90])
        elseif strcmp(pt_mode,'analysis')
            trackPlot(Xdata,[info.radius],frame)
            title('Radius [find\_cell]','FontSize',TFontSize)
        end
        if n_plots == 4; xlabel('frame number','FontSize',XFontSize); end
        
        if n_plots == 6
            subplot(sbplt(1), sbplt(2), sbplt(3) + 5)
            if strcmp(pt_mode,'data')
                if strcmp(show_fit,'ellipseDetection')
                    trackPlot(Xdata,(fits(3,:)-fits(4,:))./(fits(3,:)+fits(4,:)),frame)
                    title('Taylor Parameter [ellipse fitting]','FontSize',TFontSize)
                elseif max(strcmp(show_fit,{'unwrap','regionProps'})) == 1
                    trackPlot(Xdata,[info.uTaylorParameter],frame)
                    title('Taylor Parameter [unwrapping]','FontSize',TFontSize)
                    ylim([0 0.06])
                end
            elseif strcmp(pt_mode,'analysis')
                trackPlot(repmat(Xdata,4,1)',[info.crop]',frame)
                title('Crop box corner co-ordinates [find\_cell]','FontSize',TFontSize), ylabel('px','FontSize',YFontSize)
                legend('x0','x1','y0','y1')
            end
            %xlabel('Frame number','FontSize',XFontSize)
            
            subplot(sbplt(1), sbplt(2), sbplt(3) + 6)
            if strcmp(pt_mode,'data')
                if strcmp(show_fit,'ellipseDetection')
                    
                    %     trackPlot(Xdata, [info.MajorAxisLength]/2 - [info.radius], frame)
                    trackPlot(Xdata, fits(5,:), frame)
                    ylim([-90 90]), title('Angle of long axis above horizontal [ellipse fitting]','FontSize',TFontSize), ylabel('Angle (^o)','FontSize',YFontSize)
                elseif max(strcmp(show_fit,{'unwrap','regionProps'})) == 1
                    trackPlot(Xdata, (180/pi)*[info.uOrientation], frame, 'x')
                    ylim([-90 90])
                    title('Angle of long axis above horizontal [unwrapping]','FontSize',TFontSize), ylabel('Angle (^o)','FontSize',YFontSize)
                end
            elseif strcmp(pt_mode,'analysis')
                trackPlot(repmat(Xdata,2,1)',[info.find_fails; info.seg_fails]',frame)
                title('Success metrics [find\_cell, seg\_cell]','FontSize',TFontSize)
                legend('find','seg')
                ylim([-1 3])
            end
            xlabel('Frame number','FontSize',XFontSize)
        end
    end
    if makevid == 0 && flythrough == 0
        input('')
    elseif makevid == 0 && flythrough == 1
        pause(p_time)
    elseif makevid ==1
        drawnow
        fr = getframe(h);
        im = frame2im(fr);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if frame == 1
            imwrite(imind,cm,strcat(fname,'.gif'),'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,strcat(fname,'.gif'),'gif','WriteMode','append');
        end
    elseif makevid == 2
        drawnow
        fr = getframe(h);
        writeVideo(v, fr);
    end
end
if makevid == 2
    close(v);
    clear fr v
elseif makevid == 1
    clear im cm fr imind
end
%% These plots are a bit nicer to show people
figure(19)
subplot(221)
plot([info.TaylorParameter])
title('Taylor Parameter'), xlabel('Frame number')
ylim([0, 1.1*max([info.TaylorParameter])])
subplot(222)
plot([info.Area].*(0.07^2))
title('Area'), xlabel('frame number'), ylabel('\mu m^2')
ylim([0, 1.1*max([info.Area])]*(0.07^2))
subplot(223)
plot([info.Orientation])
title('Angle of long axis above horizontal')
xlabel('frame numer'), ylabel('Angle (degrees)')
ylim([-90 90])
subplot(224)
plot(centroids(1:2:end)*0.07,'b')
hold on
plot(centroids(2:2:end)*0.07,'k')
ylim([0 1920]*0.07), xlabel('frame number'), ylabel('\mu m')
title('Centroids of mask')
%% Show the unwrapped data for a single frame (you need to run parts of unwrap_cell manually)
for fr  =1:length(Imstack{1})
    info(fr).uMajorAxisLength = fits(1,fr);
    info(fr).uMinorAxisLength = fits(2,fr);
    info(fr).uOrientation = fits(3,fr);
end
figure(99) % I got 99 problems and a figure ain't one
imagesc(unwrapped(:,:,frs(1)))
hold on
plot(Ia(:,:,frs(1)),'b')
plot(FitEqn(info(frs(1)).uMajorAxisLength, info(frs(1)).uMinorAxisLength, ...
    info(frs(1)).uOrientation, 1:360),'k--')
hold off

%% View summary plots
% Open all plots from a given set.
Dset = 'LS174T_normoxia';
%Posn = [-3523 1100 568 434];
Posn = [-3289        1075         568         434];
plots = dir('/home/ppxwh2/Documents/data/OpTrap/processing_plots');
Offset = 10; % Files are opened in alphabetical order. Use this to shift to later files
N_max = 10;
N_plots = 0;
for item = plots'
    if N_plots < N_max
        if length(item.name) > 4
            if strcmp(item.name(end-3:end),'.fig') && contains(item.name,Dset)
                if Offset <= 0
                    %disp(item.name)
                    open(strcat('~/Documents/data/OpTrap/processing_plots/',item.name));
                    drawnow;
                    h = gcf;
                    h.Name = item.name(1:end-4);
                    Posn = Posn + [25 -25 0 0];
                    h.Position = Posn;
                    subplot(3,2,2)
                    ylim([0,0.1])
                    N_plots = N_plots + 1;
                else
                    Offset = Offset - 1;
                end
            end
        end
    end
end
%% For putting a single frame with fitting
[file, path] = uiputfile('*.png','Save your images');
frames = [50, 700];
figure(14)
for frame = frames
    aa = imagesc(Imstack{1}{frame,1});
    colormap(aa.Parent,'gray')
    hold on
    axis image off
    % Using unwrap cell fitting
    plot(info(frame).uMajorAxisLength .* cos(theta) .* cos(info(frame).uOrientation) ...
        - info(frame).uMinorAxisLength .* sin(theta) .* sin(-info(frame).uOrientation)...
        + centres(1,frame) + offset(2,frame),... % x values end here
        info(frame).uMajorAxisLength .* cos(theta) .* sin(-info(frame).uOrientation) ...
        + info(frame).uMinorAxisLength .* sin(theta) .* cos(info(frame).uOrientation)...
        + centres(2,frame) + offset(3,frame),'k--','LineWidth',2)
    plot(centres(1,frame) + offset(2,frame),centres(2,frame) + offset(3,frame),'kx')
    ax = gca;
    ax.Children(1).LineWidth = 6;
    ax.Children(1).MarkerEdgeColor = 'red';
    ax.Children(1).MarkerSize = 12;
    ax.Children(2).LineWidth = 6;
    ax.Children(2).Color = [1 0 0];
    Fr = getframe(ax);
    SaveFile = [file(1:end-4) '_fr' num2str(frame) '.png'];
    imwrite(Fr.cdata,[path SaveFile])
end
%% Function trackPlot used above
function trackPlot(xdata, ydata, idx, varargin)
    hold off
    if nargin == 3;     plot(xdata, ydata)
    elseif nargin == 4; plot(xdata, ydata, varargin{:})
    elseif nargin == 5; error('Too many nargins');
    end
    hold on
    plot(xdata(idx),ydata(idx),'rx')

end