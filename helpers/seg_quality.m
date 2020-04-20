%% Notes:
% This needs to be turned into a few functions - a loader function, a
% figure function with proper functionality (overlays, 2-6 plots, choice on
% what the plots are)
%% Load data
CellType = 'HL60';
Set = 'normoxia';
Num = '17';
global Imstack info meta;

[Imstack, info, meta] = LoadImstackInfoMeta(CellType,Set,Num);

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
n_plots = 6;                % Number of plots - 6 includes ellipse fitting results
pt_mode = 'data';           % Analysis or data or unwrap - do you want to look at the data, or analyse why it isn't working, or just show unwrapped data
frs =680:5:750;                 % Frames to display

p_time = 0.25;              % Time to pause on each frame when showing as movie
makevid = 0;                % Set to 1 to make animated gif or 2 to make avi
flythrough = 1;             % Play as movie instead of requiring user input to step through
run_unwrap = 1;             % Run unwrap cell before displaying data (refreshes content of unwrapped, fits, Ia)
UnwrapOpts = {'sc_up',1.8,'ifNaN','mean','sc_down',0.35,'UseGradient',true}; % Options for unwrap

% Font sizes for axes and titles
FontSizes.XFontSize = 12;
FontSizes.YFontSize = 12;
FontSizes.TFontSize = 12;

if isempty(whos('offset')); run_unwrap = 1; end
if run_unwrap && meta.find_cell_v; [u_fits, unwrapped, Ia, FitEqn, offset] = unwrap_cell_v2(Imstack, [info.centres] , [info.radius],UnwrapOpts{:}); 
elseif run_unwrap; [UnwrapFits, ~, FitEqn, offset, FitErrs] = unwrap_cell_v4(Imstack, [info.mCentres] , repmat(100,1,size(Imstack{1},1)),UnwrapOpts{:});end %#ok<UNRCH>

if n_plots == 4; sbplt = [4 2 4];
elseif n_plots == 6; sbplt = [6 2 6]; end

filename = strsplit(info(1).filepath,{'/','.','_'});
fname = strjoin({'seg',filename{9:15}},'_'); % Filename base for saving - seg_(original filename)
h = figure(13);

centroids = [info.Centroid];
theta = 0:0.01:2*pi;
Xdata = 1:size(info,2);
if run_unwrap; info = H_UpdateInfoUfits(info, u_fits); end
if meta.seg_cell_v == 5; fits = [info.ellipse_fits]; end
if strcmp(pt_mode,'analysis'); crops = [info.crop]; end
if strcmp(show_mask,'initial')
    [rr, cc] = meshgrid(1:size(Imstack{1}{1,1},2),...
        1:size(Imstack{1}{1,1},1));
    centres = [info.centres];
elseif strcmp(show_fit,'unwrap') && meta.line_maxima_v
    centres = [info.mCentres];
end

if makevid == 2; v = VideoWriter(strcat(fname,'.avi')); v.FrameRate = 4; open(v); end
for frame = frs
    if length(whos('unwrapped','Ia','fiteqn')) == 3
        subplot(2,2,2)
        PlotUnwrappedAndFit(unwrapped, Ia, info, FitEqn, FontSizes,frame)
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
    title(['Segmentation and fitting: frame ' num2str(frame)],'FontSize',FontSizes.TFontSize+4)
    
    % Draw fitted ellipse
    if strcmp(show_fit, 'regionProps')
        % Using regionprops fitting
        PlotEllipseOverlay(info(frame).MajorAxisLength, info(frame).MinorAxisLength, ...
            info(frame).Orientation,centroids(frame*2-1:frame*2))
    elseif strcmp(show_fit, 'ellipseDetection')
        % Using ellipseDetection fitting
        PlotEllipseOverlay(2 * info(frame).fits(3), 2 * info(frame).fits(4),...
            info(frame).fits(5),info(frame).fits(1:2))
    elseif strcmp(show_fit, 'unwrap')
        % Using unwrap cell fitting
        PlotEllipseOverlay(2*info(frame).uMajorAxisLength, 2*info(frame).uMinorAxisLength,...
            info(frame).uOrientation, centres(:,frame) + info.uOffset(2:3,frame))
    elseif strcmp(show_fit,'linemax')
        plot(Centres(1,frame),Centres(2,frame),'kx')
    end
    
    % Draw plots below
    if n_plots
        subplot(sbplt(1), sbplt(2), sbplt(3) + 1)
        cla
        N_PlotTaylorRP(Xdata, info, frame, FontSizes)
        
        subplot(sbplt(1), sbplt(2), sbplt(3) + 2)
        cla
        N_PlotAreaRP(Xdata, info, frame, FontSizes)
        
        subplot(sbplt(1), sbplt(2), sbplt(3) + 3)
        cla
        N_PlotCentreRP(centroids, frame, length(Imstack{1}{1,1}), FontSizes)
        
        if n_plots == 4; xlabel('frame number'); end
        
        subplot(sbplt(1), sbplt(2), sbplt(3) + 4)
        cla
        if strcmp(pt_mode,'data')
            N_trackPlot(Xdata,[info.Orientation],frame,'x')
            title('Angle of long axis above horizontal [regionprops]','FontSize',FontSizes.TFontSize)
            ylabel('Angle (^{o})','FontSize',FontSizes.YFontSize)
            ylim([-90 90])
        elseif strcmp(pt_mode,'analysis')
            N_trackPlot(Xdata,[info.radius],frame)
            title('Radius [find\_cell]','FontSize',FontSizes.TFontSize)
        end
        if n_plots == 4; xlabel('frame number','FontSize',FontSizes.XFontSize); end
        
        if n_plots == 6
            subplot(sbplt(1), sbplt(2), sbplt(3) + 5)
            cla
            if strcmp(pt_mode,'data')
                if strcmp(show_fit,'ellipseDetection')
                    N_trackPlot(Xdata,(fits(3,:)-fits(4,:))./(fits(3,:)+fits(4,:)),frame)
                    title('Taylor Parameter [ellipse fitting]','FontSize',FontSizes.TFontSize)
                elseif max(strcmp(show_fit,{'unwrap','regionProps'})) == 1
                    N_trackPlot(Xdata,[info.uTaylorParameter],frame)
                    title('Taylor Parameter [unwrapping]','FontSize',FontSizes.TFontSize)
                    ylim([0 0.06])
                end
            elseif strcmp(pt_mode,'analysis')
                N_trackPlot(repmat(Xdata,4,1)',[info.crop]',frame)
                title('Crop box corner co-ordinates [find\_cell]','FontSize',FontSizes.TFontSize), ylabel('px','FontSize',FontSizes.YFontSize)
                legend('x0','x1','y0','y1')
            end
            %xlabel('Frame number','FontSize',FontSizes.XFontSize)
            
            subplot(sbplt(1), sbplt(2), sbplt(3) + 6)
            if strcmp(pt_mode,'data')
                if strcmp(show_fit,'ellipseDetection')
                    
                    %     trackPlot(Xdata, [info.MajorAxisLength]/2 - [info.radius], frame)
                    N_trackPlot(Xdata, fits(5,:), frame)
                    ylim([-90 90]), title('Angle of long axis above horizontal [ellipse fitting]','FontSize',FontSizes.TFontSize), ylabel('Angle (^o)','FontSize',FontSizes.YFontSize)
                elseif max(strcmp(show_fit,{'unwrap','regionProps'})) == 1
                    N_trackPlot(Xdata, (180/pi)*[info.uOrientation], frame, 'x')
                    ylim([-90 90])
                    title('Angle of long axis above horizontal [unwrapping]','FontSize',FontSizes.TFontSize), ylabel('Angle (^o)','FontSize',FontSizes.YFontSize)
                end
            elseif strcmp(pt_mode,'analysis')
                N_trackPlot(repmat(Xdata,2,1)',[info.find_fails; info.seg_fails]',frame)
                title('Success metrics [find\_cell, seg\_cell]','FontSize',FontSizes.TFontSize)
                legend('find','seg')
                ylim([-1 3])
            end
            xlabel('Frame number','FontSize',FontSizes.XFontSize)
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
    info(fr).uTaylorParameter = (fits(1,fr)-fits(2,fr))/(fits(1,fr)+fits(2,fr));
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
Dset = 'HL60_with_drugs';
RePosn = false;

%Posn = [-3523 1100 568 434];
Posn = [-3289        1075         568         434];
plots = dir('~/Documents/data/OpTrap/processing_plots');
Offset = 0; % Files are opened in alphabetical order. Use this to shift to later files
N_max = 20;
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
                    if RePosn
                        Posn = Posn + [25 -25 0 0];
                        h.Position = Posn;
                    end
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
%% Functions used above
function N_PlotTaylorRP(Xdata, info, frame, FontSizes)
N_trackPlot(Xdata,[info.TaylorParameter],frame)
title('Taylor Parameter [regionprops]','FontSize',FontSizes.TFontSize)
YMax = 1.1*max([info.TaylorParameter]);
if isnan(YMax); YMax = 1; end
ylim([0, YMax])
end

function N_PlotAreaRP(Xdata, info, frame, FontSizes)
N_trackPlot(Xdata,(0.07^2)*[info.Area],frame)
title('Area [regionprops]','FontSize',FontSizes.TFontSize), ylabel('\mu m^2','FontSize',FontSizes.YFontSize)
YMax = 1.1*max([info.Area])*(0.07^2);
if isnan(YMax); YMax = 1; end
ylim([0, YMax])
end
        
function N_PlotCentreRP(centroids, frame, FrSize, FontSizes)
plot(centroids(1:2:end)*0.07)
hold on
plot(centroids(2:2:end)*0.07)
plot(frame,centroids(2*frame - 1)*0.07,'rx')
plot(frame,centroids(2*frame)*0.07, 'rx')
ylim([0 FrSize]*0.07), ylabel('\mu m','FontSize',FontSizes.YFontSize)
title('Centroids of mask [regionprops]','FontSize',FontSizes.TFontSize)
end



function N_trackPlot(xdata, ydata, idx, varargin)
    hold off
    if nargin == 3;     plot(xdata, ydata)
    elseif nargin == 4; plot(xdata, ydata, varargin{:})
    elseif nargin == 5; error('Too many nargins');
    end
    hold on
    plot(xdata(idx),ydata(idx),'rx')

end