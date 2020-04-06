function SubplotTitle(Imstack, info, Frs, Titles, FSizes)
% Needs some work - should take plot commands as a cell array and then use
% evalin to do plotting

Frames = [Frs; Frs];
% M, N, and n should be calculated from a given N_plots
M = 11;
N = 11;
subplot(M,N,1:N)
ax = gca;
ax.Color = ax.Parent.Color;
ax.XColor = ax.Parent.Color;
ax.YColor = ax.Parent.Color;
title([CellType ' ' Set ' '],'FontSize',FSizes.Ttl1)

for n = 0:3
    % The workhorse of the loop - the list of subplots covered
    V = floor(n/2) * N*(M-1)/2 + mod(n,2) * (N+1)/2 + (1:(N-1)/2) + (2*N:N:N*(-3+M)/2)';
    subplot(M,N,reshape(V,1,[]))
    imagesc(Imstack{1}{Frames(n+1),1})
    title(Titles{n+1},'FontSize',FSizes.Ttl2)
    axis image off, hold on
    
    
    PlotEllipseOverlay(2 * info(fr).uMajorAxisLength, 2*info(fr).uMinorAxisLength,...
        info(fr).uOrientation, info(fr).mCentres + info(fr).uOffset(2:3))
    
end 

end