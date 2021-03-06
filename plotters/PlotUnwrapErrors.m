function PlotUnwrapErrors(CellType, Set, Num, varargin)
%% [Imstack, info, meta] = plot_unwrap_errors(CellType, Set, Num, varargin)
%% Plot deformation for one dataset with errors from std of relaxed deformation
% This needs a proper input parser soon

global Imstack info meta

if nargin > 3
    [SvFig, SvPng] = varargin{:};
else
    SvFig = false;
    SvPng = false;
end

% Check correct info, meta and Imstack are loaded
LoadImstackInfoMeta(CellType,Set,Num);

%% Plotting variables
FSize = 16;
Frs = 1:90;

Tdata = linspace(0,meta.N_Frames/100,meta.N_Frames);
%% Calculate error from standard deviation of deformation from select frames
DErrs = repmat(std([info(Frs).uTaylorParameter],0,2),1,meta.N_Frames);
%figure(Fh)
%clf
hold on
errorbar(Tdata, [info.uTaylorParameter],DErrs,'.','MarkerEdgeColor','r');
% ylim([0 0.06])
xlabel('Time (s)','FontSize',FSize)
ylabel('Deformation (Taylor Paramater: scale 0 to 1)','FontSize',FSize)
title({'Deformation with effors from standard' 'deviation of relaxed cell deformation' ['Set: ' strjoin({CellType,Set,Num})]},...
    'FontSize',FSize)
if SvFig
    saveas(gcf,['~/fig/EQ/D_' strjoin({CellType, Set, Num},'_') '_std_D.fig'])
end
if SvPng
    saveas(gcf,['~/png/EQ/D_' strjoin({CellType, Set, Num},'_') '_std_D.png'])
end

end