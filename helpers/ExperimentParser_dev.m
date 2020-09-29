%% Experiment parser
% Use position from find_cell to determine when the cell is being dragged
% (and in which direction?)

CellType = 'HL60';
Set = 'normoxia';
Num = '19';

global Imstack info meta
LoadImstackInfoMeta(CellType, Set, Num, false, '', 'find5',true)

%%


% Crop start and end to remove ringing
CropFrames = 50;
% How many frames to use in calculating standard deviation
STDFrames = 800;
% Threshold scaling for block ends (higher is less sensitive)
ThreshFactor = 1;
% Minimum length for a block (below this is discarded for being noise)
MinBlockSize = 100;
% Min and Max number of blocks
MinBlocks = 3;
MaxBlocks = 5;
% Maximum number of iterations
MaxRunCount = 1000;

% Get centres array and low-pass filter it
Centres = [info.centres]';
LPCentres = lowpass(Centres,1,100);
% Get number of frames from centres array
N_Frames = size(Centres,1);
% Elementwise difference (D(i) = x(i) - x(i-1)) thresholded relative to
% standard deviation
DiffCentres = diff(LPCentres(CropFrames:end-CropFrames,:),1,1);
BlocksX = [0;0];
BlocksY = [0;0];
BlockLengthsX = 0;
BlockLengthsY = 0;
RunCount = 0;

while (length(BlocksY) < MinBlocks || length(BlocksY) > MaxBlocks || min(BlockLengthsY) < MinBlockSize) && RunCount < MaxRunCount
    if RunCount == 0
    elseif length(BlocksY) < MinBlocks && ThreshFactor > 1
        ThreshFactor = 0.9 * ThreshFactor;
    elseif length(BlocksY) > MaxBlocks && ThreshFactor < 4
        ThreshFactor = 1.1 * ThreshFactor;
    elseif length(BlocksY) < MinBlocks
        ThreshFactor = 2;
        STDFrames = ceil(STDFrames * 1.2);
    elseif length(BlocksY) >= MaxBlocks
        ThreshFactor = 2;
        STDFrames = ceil(STDFrames * 0.8);
    elseif sum(BlockLengthsY < MinBlockSize) <= 2
        NumSmallBlocks = sum(BlockLengthsY < MinBlockSize);
        warning(['Cannot prevent ' num2str(NumSmallBlocks) ' small blocks, breaking'])
        break
    else
        error('The situation you said wouldn''t happen... happened. Somehow.')
    end
    RunCount = RunCount + 1;
    
    STDDiffCentres = std(DiffCentres(1:STDFrames,:));
    ThresholdedDiff = abs(DiffCentres) > STDDiffCentres .* ThreshFactor;
    % If D(i+1) > D(i) you're in a new block of changing position
    ShiftedR = zeros(size(ThresholdedDiff),'logical');
    ShiftedR(2:end,:) = ThresholdedDiff(1:end-1,:);
    BlockStart = ShiftedR > ThresholdedDiff;
    % If D(i-1) > D(i) you're in a new block of equilibrium
    ShiftedL = zeros(size(ThresholdedDiff),'logical');
    ShiftedL(1:end-1,:) = ThresholdedDiff(2:end,:);
    BlockEnd = ShiftedL > ThresholdedDiff;
    % First block starts at 1, the rest are shifted by being cropped
    BlocksStartX = [1, find(BlockStart(:,1))' + CropFrames];
    BlocksStartY = [1, find(BlockStart(:,2))' + CropFrames];
    % Last block ends at N_Frames
    BlocksEndX = [find(BlockEnd(:,1))' + CropFrames, meta.N_Frames];
    BlocksEndY = [find(BlockEnd(:,2))' + CropFrames, meta.N_Frames];
    % Blocks ready for output
%     try
        BlocksX = [BlocksStartX; BlocksEndX];
        BlocksY = [BlocksStartY; BlocksEndY];
        BlockLengthsX = diff(BlocksX);
        BlockLengthsY = diff(BlocksY);
%     catch
%         BlocksX = [];
%         BlocksY = [];
%         BlockLengthsX = [];
%         BlockLengthsY = [];
%     end
    
    

end

if isempty(BlocksX)
    error('Empty blocks array returned')
end
BlocksX = BlocksX(:, BlockLengthsX > MinBlockSize);
BlocksY = BlocksY(:, BlockLengthsY > MinBlockSize);
meta.Experiment.BlocksX = {BlocksX};
meta.Experiment.BlocksY = {BlocksY};
% %%
% Starts = repmat(BlocksStartX', 1, length(BlocksEndX));
% Ends = repmat(BlocksEndX, length(BlocksStartX), 1);
% X = Starts == Ends;
% sum(X,'all')
% [I, J] = find(X);
% BlocksStartX(I)
% BlocksEndX(J)

%%
clf
subplot(2,1,1)
N_PlotBlocks(Centres, meta.N_Frames, CropFrames, BlocksX, BlocksY, STDFrames + CropFrames) %

subplot(2,1,2)
plot(CropFrames+(1:length(DiffCentres)),abs(DiffCentres).*[1, -1])
hold on
plot([0 0; length(DiffCentres).* [1 1] ],[1 -1; 1 -1] .* STDDiffCentres .* ThreshFactor)
% subplot(3,1,3)
% plot(BlockStart,'x')

function N_PlotBlocks(Centres, N_Frames, Xshift, BlocksX, BlocksY, EndOfSTD)
% Plots centres and blocks on current axis
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
plot([1 1].* EndOfSTD, Y(EndOfSTD,:)','k.')

legend('X centre','Y centre','Last frame used in \sigma calculation')
end
