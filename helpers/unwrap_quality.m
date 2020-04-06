%% Compare unwrap fits
% Show fit on cell and fit from different frame on cell


% Load
CellType = 'HL60';
Set = 'normoxia';
Num = '18';

args = N_TidyLoader(CellType, Set, Num);
if length(args)>1
    [Imstack, info, meta, Ia, FitEqn] = args{:};
end

Frs = [1, 900];

fh = figure(37);
colormap gray
for fr = Frs
    subplot(2,length(Frs), find(Frs==fr))
    imagesc(Imstack{1}{fr,1})
    axis image off, hold on
    PlotEllipseOverlay(2 * info(fr).uMajorAxisLength, 2*info(fr).uMinorAxisLength,...
        info(fr).uOrientation, info(fr).mCentres + info(fr).uOffset(2:3))
    
    subplot(2, length(Frs), length(Frs) + find(Frs==fr))
    imagesc(Imstack{1}{fr,1})
    axis image off, hold on
    PlotEllipseOverlay(2 * info(fr).uMajorAxisLength, 2*info(fr).uMinorAxisLength,...
        info(fr).uOrientation, info(fr).mCentres + info(fr).uOffset(2:3))
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
end
end