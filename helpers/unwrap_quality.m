%% Compare unwrap fits
% Show fit on cell and fit from different frame on cell

% Load
CellType = 'HL60';
Set = 'normoxia';
Num = '18';

% Display options
Frs = [1, 500]; % Which frames

% Fontsizes
FSizes.Ttl1 = 16; % Titles
FSizes.Ttl2 = 14;
FSizes.XL1 = 14; % YLabels
FSizes.XL2 = 10;
FSizes.YL1 = 14; % XLabels
FSizes.YL2 = 10;

args = N_TidyLoader(CellType, Set, Num);
if length(args)>1
    [Imstack, info, meta, Ia, FitEqn] = args{:};
end

%% Show overlaid images and fits
fh = figure(37);
colormap gray
for fr = Frs
    subplot(2,length(Frs), find(Frs==fr))
    imagesc(Imstack{1}{fr,1})
    axis image off, hold on
    PlotEllipseOverlay(2 * info(fr).uMajorAxisLength, 2*info(fr).uMinorAxisLength,...
        info(fr).uOrientation, info(fr).mCentres + info(fr).uOffset(2:3))
    title(['Frame ' num2str(fr) ' with fit'])
    
    subplot(2, length(Frs), length(Frs) + find(Frs==fr))
    imagesc(Imstack{1}{fr,1})
    axis image off, hold on
    PlotEllipseOverlay(2 * info(Frs(Frs~=fr)).uMajorAxisLength, 2*info(Frs(Frs~=fr)).uMinorAxisLength,...
        info(Frs(Frs~=fr)).uOrientation, info(Frs(Frs~=fr)).mCentres + info(Frs(Frs~=fr)).uOffset(2:3))
    title(['Frame ' num2str(fr) ' with fit from frame ' num2str(Frs(Frs~=fr))])
end

%% Same but with nice title


Titles = {['Frame ' num2str(Frs(1)) ' with fit from frame ' num2str(Frs(2))]; 
    ['Frame ' num2str(Frs(2)) ' with fit from frame ' num2str(Frs(2))]; 
    ['Frame ' num2str(Frs(1)) ' with fit from frame ' num2str(Frs(1))]; 
    ['Frame ' num2str(Frs(2)) ' with fit from frame ' num2str(Frs(1))]};
Frames = [Frs; Frs];
M = 11;
N = 11;

fh = figure(38);
clf
subplot(M,N,1:N)
ax = gca;
ax.Color = ax.Parent.Color;
ax.XColor = ax.Parent.Color;
ax.YColor = ax.Parent.Color;
title([CellType ' ' Set ' example fits with CI'],'FontSize',FSizes.Ttl1)

for n = 0:3
    % The workhorse of the loop - the list of subplots covered
    V = floor(n/2) * N*(M-1)/2 + mod(n,2) * (N+1)/2 + (1:(N-1)/2) + (2*N:N:N*(-3+M)/2)';
    subplot(M,N,reshape(V,1,[]))
    imagesc(Imstack{1}{Frames(n+1),1})
    
    PlotEllipseOverlay(2 * info(fr).uMajorAxisLength, 2*info(fr).uMinorAxisLength,...
        info(fr).uOrientation, info(fr).mCentres + info(fr).uOffset(2:3))
    
    title(Titles{n+1},'FontSize',FSizes.Ttl2)
    axis image off, hold on
    xlabel('
end 


%%
function [varargout] = N_TidyLoader(CellType, Set, Num)
try
    evalin('base','compare_info_meta_imstack(info, meta, Imstack)')
    varargout = {true};
catch
    [Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num);
    UnwrapOpts = {'UseGradient',true};
    if ~meta.line_maxima_v
        [u_fits, ~, Ia, FitEqn, ~] = unwrap_cell_v2(...
            Imstack, [info.centres] , [info.radius],UnwrapOpts{:});
    else
        [u_fits, ~, Ia, FitEqn, ~, ~] = unwrap_cell_v2(...
            Imstack, [info.mCentres] , repmat(100,1,size(Imstack{1},1)),UnwrapOpts{:});
    end
    info = H_UpdateInfoUfits(info, u_fits);
    varargout = {Imstack, info, meta, Ia, FitEqn};
    %}
end
end