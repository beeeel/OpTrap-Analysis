%% Compare unwrap fits
% Show fit on cell and fit from different frame on cell

% Load
CellType = 'MV411';
Set = 'normoxia';
Num = '1';

% Display options
Frs = [1, 900]; % Which frames

% Fontsizes
FSizes.Ttl1 = 20; % Titles
FSizes.Ttl2 = 14;
FSizes.XL1 = 14; % YLabels
FSizes.XL2 = 10;
FSizes.YL1 = 14; % XLabels
FSizes.YL2 = 10;

global Imstack info meta
LoadImstackInfoMeta(CellType, Set, Num);

%% Show overlaid images and fits with confidence interval
% Fits are rotated and translated to orientation and centre of cell.
% Confidence interval from standard deviation of initial deformation.
SaveFig = false;
SavePng = false;

Ds = [info.uTaylorParameter];
DErrs = repmat(std(Ds(1:100)),1,length(info));
as = [info.uMajorAxisLength];
bs = [info.uMinorAxisLength];

aErrs = sqrt(2 * (as.^2 .* (DErrs./Ds).^2 .* (1 - Ds) ./ (as + bs).^2));
bErrs = sqrt(2 * (bs.^2 .* (DErrs./Ds).^2 .* (1 + 3 * Ds) ./ (as + bs).^2));
FSpec = 'D = %.4f ± %.4f';
Titles = {{['Frame ' num2str(Frs(1)) ' with fit from frame ' num2str(Frs(1))], sprintf(FSpec ,info(Frs(1)).uTaylorParameter, DErrs(1))}; 
    {['Frame ' num2str(Frs(1)) ' with fit from frame ' num2str(Frs(2))], sprintf(FSpec ,info(Frs(2)).uTaylorParameter, DErrs(1))}; 
    {['Frame ' num2str(Frs(2)) ' with fit from frame ' num2str(Frs(1))], sprintf(FSpec ,info(Frs(1)).uTaylorParameter, DErrs(1))}; 
    {['Frame ' num2str(Frs(2)) ' with fit from frame ' num2str(Frs(2))], sprintf(FSpec ,info(Frs(2)).uTaylorParameter, DErrs(1))}};
Frames = [Frs; Frs];
M = 11;
N = 11;

fh = figure(38);
clf
colormap gray
subplot(M,N,1:N)
ax = gca;
ax.Color = ax.Parent.Color;
ax.XColor = ax.Parent.Color;
ax.YColor = ax.Parent.Color;
title({[CellType ' ' Set ' ' Num ' example fits with CI'], ['Showing fitted axis + confidence interval. ' num2str(meta.N_Frames) ' frames.']},...with errors from','standard deviation of initial deformation'},...
    'FontSize',FSizes.Ttl1)

for n = 0:3
    fr2 = Frs(mod(n,2)+1); 
    fr1 = Frs(ceil((n+1)/2));
    % The workhorse of the loop - the list of subplots covered
    V = floor(n/2) * N*(M-1)/2 + mod(n,2) * (N+1)/2 + (1:(N-1)/2) + (2*N:N:N*(-3+M)/2)';
    subplot(M,N,reshape(V,1,[]))
    imagesc(Imstack{1}{fr1,1})
    axis image off, hold on
    % The offset seemed to make it fit wrong. This has been "hacked out".
    PlotEllipseOverlay(2 * info(fr2).uMajorAxisLength, ...
        2*info(fr2).uMinorAxisLength,...
        info(fr1).uOrientation, info(fr1).mCentres + 0* info(fr1).uOffset(2:3),'k','LineWidth',3)
    PlotEllipseOverlay(2 * (info(fr2).uMajorAxisLength + aErrs(fr2)), ...
        2 * (info(fr2).uMinorAxisLength - bErrs(fr2) ),...
        info(fr1).uOrientation, info(fr1).mCentres + 0* info(fr1).uOffset(2:3),'r:','LineWidth',2)
    PlotEllipseOverlay(2 * (info(fr2).uMajorAxisLength - aErrs(fr2)), ...
        2 * (info(fr2).uMinorAxisLength + bErrs(fr2) ),...
        info(fr1).uOrientation, info(fr1).mCentres + 0* info(fr1).uOffset(2:3),'b--','LineWidth',2)
    title(Titles{n+1},'FontSize',FSizes.Ttl2)
end

SaveFigPng(['fits_CI_std_D_' strjoin({CellType,Set,Num},'_')],'EQ',SaveFig,SavePng)