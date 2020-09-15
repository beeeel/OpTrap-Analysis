%% Experiment parser
% Use position from find_cell to determine when the cell is being dragged
% (and in which direction?)

CellType = 'HL60';
Set = 'normoxia';
Num = '7';

global Imstack info meta
LoadImstackInfoMeta(CellType, Set, Num, false, '', 'find5',true)

%%
% Crop start and end to remove ringing
CropFrames = 50;
% Threshold scaling for block ends (higher is less sensitive)
ThreshFactor = 2;
% Get centres array and low-pass filter it
Centres = [info.centres]';
LPCentres = lowpass([info.centres]',1,100);
% Elementwise difference (D(i) = x(i) - x(i-1)) thresholded relative to
% standard deviation
DiffCentres = diff(LPCentres(CropFrames:end-CropFrames,:),1,1);
STDDiffCentres = std(DiffCentres);
ThresholdedDiff = abs(DiffCentres) > STDDiffCentres .* ThreshFactor;
% If D(i+1) > D(i) you're in a new block of changing position
ShiftedR = zeros(size(ThresholdedDiff),'logical');
ShiftedR(2:end) = ThresholdedDiff(1:end-1);
BlockStart = ShiftedR > ThresholdedDiff;
% If D(i-1) > D(i) you're in a new block of equilibrium
ShiftedL = zeros(size(ThresholdedDiff),'logical');
ShiftedL(1:end-1) = ThresholdedDiff(2:end);
BlockEnd = ShiftedL > ThresholdedDiff;
% First block starts at 1, the rest are shifted by being cropped
BlocksStartX = [1, find(BlockStart(:,1))' + CropFrames];
BlocksStartY = [1, find(BlockStart(:,2))' + CropFrames];
% Last block ends at N_Frames
BlocksEndX = [find(BlockEnd(:,1))' + CropFrames, meta.N_Frames];
BlocksEndY = [find(BlockEnd(:,2))' + CropFrames, meta.N_Frames];
% Blocks ready for output
BlocksX = [BlocksStartX; BlocksEndX];
BlocksY = [BlocksStartY; BlocksEndY];
BlockLengthsX = diff(BlocksX);
BlockLengthsY = diff(BlocksY);
BlocksX = BlocksX(:, BlockLengthsX > 100);
BlocksY = BlocksY(:, BlockLengthsY > 100);

info.Experiment.BlocksX = {BlocksX};
info.Experiment.BlocksY = {BlocksY};
%%
clf
subplot(2,1,1)
N_PlotBlocks(Centres, meta.N_Frames, CropFrames, BlocksX, BlocksY) %

subplot(2,1,2)
plot(CropFrames+(1:length(DiffCentres)),abs(DiffCentres).*[1, -1])
hold on
plot([0 0; length(DiffCentres).* [1 1] ],[1 -1; 1 -1] .* STDDiffCentres .* ThreshFactor)
% subplot(3,1,3)
% plot(BlockStart,'x')

function N_PlotBlocks(Centres, N_Frames, Xshift, BlocksX, BlocksY)
%% Plots centres and blocks on current axis
hold on
YShift = 100;
% Plot the blocks data
Blocks = [BlocksX BlocksY];
for i = 1:length(Blocks)
    X = Blocks(1,i);
    H = 100;
    Y = (i > size(BlocksX,2)) * (-YShift - H) + YShift/2;
    W = diff(Blocks(:,i));
    rectangle('Position',[X Y W H],'FaceColor',[.5 .5 .5],'LineWidth',1,'EdgeColor',[1 1 1])
end
% Plot the centres data
X = Xshift : N_Frames - Xshift;
Y = normalize(Centres(Xshift:end-Xshift,:),'center') + [1 -1] .* YShift;
plot(X, Y(:,1), X, Y(:,2))
end
