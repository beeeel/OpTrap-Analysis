function [Imstack, info, meta] = PlotUnwrapErrors(CellType, Set, Num, varargin)
%% [Imstack, info, meta] = plot_unwrap_errors(CellType, Set, Num, varargin)
%% Plot deformation for one dataset with errors from std of relaxed deformation
% This needs a proper input parser soon

if nargin > 3
    [SvFig, SvPng] = varargin{:};
else
    SvFig = true;
    SvPng = true;
end

[Imstack, info, meta] = LoadImstackInfoMeta(CellType,Set,Num);

% Check correct info, meta and Imstack are loaded

if ~isempty(whos('info')) && ~isempty(whos('meta')) && ~isempty(whos('Imstack'))
    compare_info_meta_imstack(info, meta, Imstack)
end

N_frames = size(Imstack{1},1);

%% Plotting variables
FSize = 16;
Frs = 1:90;

Tdata = linspace(0,9.99,size(Imstack{1},1));
%% Calculate error from standard deviation of deformation from select frames
DErrs = repmat(std([info(Frs).uTaylorParameter],0,2),1,N_frames);
%figure(Fh)
%clf
hold on
errorbar(Tdata, [info.uTaylorParameter],DErrs,'.','MarkerEdgeColor','r');
ylim([0 0.06])
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