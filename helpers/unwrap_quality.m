%% Compare unwrap fits
% Show fit on cell and fit from different frame on cell


% Load
CellType = 'HL60';
Set = 'normoxia';
Num = '17';
[Imstack, info, meta, Ia, FitEqn] = N_TidyLoader(CellType, Set, Num);

Frs = [1, 900];

fh = figure(37);
colormap gray
for fr = Frs
    subplot(2,length(Frs), find(Frs==fr))
    imagesc(Imstack{1}{fr,1})
    axis image off
    PlotEllipseOverlay(2 * info(frame).uMajorAxisLength, 2*info(frame.uMinorAxisLength),...
        info(frame).uOrientation, centres(:, frame) + offset(2:3,frame))
end


%% 
function [Imstack, info, meta, Ia, FitEqn] = N_TidyLoader(CellType, Set, Num)
try
    compare_info_meta_imstack(info, meta, Imstack)
catch
    [Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num);
    UnwrapOpts = {'UseGradient',true};
    if meta.find_cell_v
        [u_fits, ~, Ia, FitEqn, ~] = unwrap_cell_v2(...
            Imstack, [info.centres] , [info.radius],UnwrapOpts{:});
    else
        [u_fits, ~, Ia, FitEqn, ~, ~] = unwrap_cell_v2(...
            Imstack, [info.mCentres] , repmat(100,1,size(Imstack{1},1)),UnwrapOpts{:});
    end
    info = H_UpdateInfoUfits(info, u_fits);
end
end