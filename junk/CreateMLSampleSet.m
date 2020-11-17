 %% Take some images out to create training dataset
% Output images directory
SaveDir = '~/AI_Optrap/Samples/Set2/';
% Dummy class - images in this class are all 0s. 
% 'stretched' or 'relaxed'. Anything else will be ignored/treated as none
DummyClass = 'none';
% Output size expression to be evaluated within the deepest loop 
% OutSizeExp = ''; 

%Augmentation:
% Random shift - N px in random direction, applied to every image (doubles
% size of dataset)
RandomShift = 20;
% Stretch - Use interpolation to stretch fractionally. [1 2] would double
% the width, for example. Only applied to "stretched" image class.
StretchRatio = [1 1.1];

% Which datasets to take images from
Cells = {'LS174T', 'HL60', 'MV411'};
DSets = {{'normoxia','hypoxia'},{'normoxia','with_drugs'},{'normoxia','with_drugs'}};
Nums = {{1:20, 1:20},{1:30, 1:20},{1:30, 1:20}};
% How many images to take from each dataset
FramesPerSet = 150;

% Calculate how many images there will be
AllNums = [Nums{:}];
NumSets = length([AllNums{:}]);
NumIms = FramesPerSet * NumSets; % Number of unique images
TotalIms = NumIms * (1 + (RandomShift ~= 0) + 0.5 * max(StretchRatio ~= 1)); % Number of images after augmentation

SelectedIms = cell(TotalIms, 2);
Count = 0;
global Imstack info meta
StartTime = tic;

for CTidx = 1:length(Cells)
    CellType = Cells{CTidx};
    disp(['Started ' CellType])
    for DSidx = 1:length(DSets{CTidx})
        DSet = DSets{CTidx}{DSidx};
        disp(['Started ' DSet])
        for Num = Nums{CTidx}{DSidx}
            NumStr = num2str(Num);
            disp(['Loading ' NumStr])
            LoadImstackInfoMeta(CellType, DSet, NumStr, true)
            N_Frames = size(Imstack{1},1);
            %             BlockSize = meta.N_Frames./FramesPerSet;
            %             Offset = randi(BlockSize);
            %             for frame = 0:FramesPerSet-1
            for frame = 1:FramesPerSet
                Count = Count + 1;
                if frame <= FramesPerSet/2
                    Offset = 0;
                elseif N_Frames == 1000
                    Offset = 750 - FramesPerSet/2;
                elseif N_Frames == 2000
                    Offset = 1650 - FramesPerSet/2;
                else
                    error([N_Frames ' frames in this set: ' meta.filepath])
                end
                
                ImInfo = struct;
                ImInfo.CellType = CellType;
                ImInfo.Set = DSet;
                ImInfo.Num = NumStr;
                ImInfo.FrNum = frame + Offset; %frame * BlockSize + Offset;
                ImInfo.ImW = size(Imstack{1}{ImInfo.FrNum},2);
                ImInfo.ImH = size(Imstack{1}{ImInfo.FrNum},1);
                ImInfo.Augment = '';
                
                % Round up output size to nearest power of 2
                %OutSize = repmat(2^ceil(log2(max(size(Imstack{1}{1,1})))),1,2);
                % Output size = image size
%                 OutSize = size(Imstack{1}{1,1});
                % Output size = minimum of image size (square output)
                OutSize = repmat(min(size(Imstack{1}{1,1})),1,2) - RandomShift;
                % Output size = 512 * 512
                %OutSize = min([512, 512],[ImInfo.ImH, ImInfo.ImW]);
                % Output size = 256 * 256
%                 OutSize = [256, 256];
                
                SelectedIms{Count,1} = CropAndShift(...
                    Imstack{1}{ImInfo.FrNum}, OutSize, 0);
                SelectedIms{Count,2} = ImInfo;
                % If RandomShift is a variable with a value
                if exist('RandomShift','var') == 1 && RandomShift
                    Count = Count + 1;
                    SelectedIms{Count,1} = CropAndShift(...
                        Imstack{1}{ImInfo.FrNum}, OutSize, RandomShift);
                    SelectedIms{Count,2} = ImInfo;
                    SelectedIms{Count,2}.Augment = 'Translation';
                end
                % If StretchRatio is a variable with at least one non-one
                % value AND we're selecting a stretched frame
                if exist('StretchRatio','var') == 1 && Offset ~= 0 && max(StretchRatio ~= 1)
                    Count = Count + 1;
                    SelectedIms{Count,1} = CropAndStretch(...
                        Imstack{1}{ImInfo.FrNum},OutSize, StretchRatio);
                    SelectedIms{Count,2} = ImInfo;
                    SelectedIms{Count,2}.Augment = 'Stretch';
                end
            end
        end
    end
end
toc(StartTime)
Vars = whos;
fprintf('Mem size of SelectedIms cell: %g GB\n',Vars(strcmp({Vars.name},'SelectedIms')).bytes./1e9)
%% Augment - WIP
% Take a circular mask, interpolate the sample to create deformation, stick
% back on top of the image
% for Imidx = 1:NumIms
%     
% end
%% Save these files as .bmp (bitmap)
LastNum = '';
RelaxedCount = 1;
StretchedCount = 1;

% Catch some annoying errors around non-existent directories
if exist([SaveDir 'relaxed'],'dir') == 0 || exist([SaveDir 'stretched'],'dir') == 0
    system(['mkdir -p ' SaveDir 'relaxed']);
    system(['mkdir -p ' SaveDir 'stretched']);
elseif exist(SaveDir,'dir') ~= 7
    error('SaveDir is expected to be a directory')
elseif ~strcmp(SaveDir(end),'/')
    SaveDir = [SaveDir '/'];
end

for Imidx = 1:TotalIms
    II = SelectedIms{Imidx,2};
    if II.FrNum > 500 || strcmp(II.Augment,'Stretch')
        Dir = [SaveDir 'stretched/'];
        ThisImClass = 'stretched';
        % These lines (and 2 in the else) are for file names as 1,2,3...
%         FName = num2str(StretchedCount); 
%         StretchedCount = StretchedCount + 1;
    else
        Dir = [SaveDir 'relaxed/'];
        ThisImClass = 'relaxed';
%         FName = num2str(RelaxedCount);
%         RelaxedCount = RelaxedCount + 1;
    end
    if isempty(II.Augment)
        Aug = '';
    else
        Aug = II.Augment(1);
    end
    FName = strjoin({II.CellType, II.Set, II.Num, num2str(II.FrNum), Aug},'_');
%     if strcmp(DummyClass,ThisImClass)
    if true%(strcmp(ThisImClass,'relaxed') && isempty(II.Augment)) || strcmp(II.Augment,'Stretch')
        imwrite(uint8(SelectedIms{Imidx,1}),[Dir FName '.png'])
    end
%     else
%         imwrite(zeros(II.ImH,II.ImW,'uint8'),[Dir FName '.png'])
%     end
%     [FID, msg] = fopen([Dir FName '.txt'],'w');
%     fprintf(FID,'CellType\t%s\tSet\t%s\tSetNum\t%s\tFrame\t%i\tRadius\t%g\tCentre(x)\t%g\tCentre(y)\t%g',II.CellType,II.Set,II.Num,II.FrNum,II.Radius,II.Centre(1),II.Centre(2));
%     fclose(FID);
    if ~strcmp(II.Num,LastNum)
        disp(['Started files for ' II.CellType ' ' II.Set ' ' II.Num])
    end
    LastNum = II.Num;
end

function ImOut = CropAndShift(ImIn, OutSize, ShiftMag)
    %% Crop and potentially shift (in x,y) an input image, returning array of matching type
    % Get image size and calculate range of indices to take
    [ImH, ImW] = size(ImIn);
    if max(OutSize > [ImH, ImW])
        fprintf('Output size: [%i, %i].\nInput size: [%i, %i].\n',...
            OutSize(1), OutSize(2), ImH, ImW)
        error('Requested output size larger than input size')
    elseif sum(([ImH, ImW] - OutSize) < ShiftMag) > 0
        warning('Shift requested is larger than difference of sizes. Actual shift may be less than expected')
    end
    Xidx = [floor((ImW-OutSize(2))/2)+1, floor((ImW+OutSize(2))/2)];
    Yidx = [floor((ImH-OutSize(1))/2)+1, floor((ImH+OutSize(1))/2)];
    % If a shift has been requested, calculate a shift in a random
    % direction
    if abs(ShiftMag) > 0
        Theta = 2*pi*rand(1);
        Xidx = Xidx + floor(ShiftMag * cos(Theta));
        Yidx = Yidx + floor(ShiftMag * sin(Theta));
        % Shift to ensure first index is within source image
        Xidx = Xidx + (Xidx(1) < 1) * (1 - Xidx(1));
        Yidx = Yidx + (Yidx(1) < 1) * (1 - Yidx(1));
        % Shift to ensure last index is within source image
        Xidx = Xidx + (Xidx(2) > ImW) * (ImW - Xidx(2));
        Yidx = Yidx + (Yidx(2) > ImH) * (ImH - Yidx(2));
    end
    if Xidx(2) > ImW
        disp('a');
    end
    ImOut = ImIn(Yidx(1):Yidx(2),Xidx(1):Xidx(2));
end

function ImOut = CropAndStretch(ImIn, OutSize, StretchRatio)
    %% Crop and stretch an input image, returning array of the same type
    [ImH, ImW] = size(ImIn);
    if max(OutSize > [ImH, ImW])
        fprintf('Output size: [%i, %i].\nInput size: [%i, %i].\n',...
            OutSize(1), OutSize(2), ImH, ImW)
        error('Requested output size larger than input size')
    end
    % Input image space
    X = 1:ImW;
    Y = 1:ImH;
    % Sampling points - start at the centre (ImW/2 e.g.) subtract half the
    % stretched space (which is smaller than image space, hence divide) to
    % find the input space co-ordinate for the edge of the output image.
    Xq = linspace(floor(ImW/2 - OutSize(2)/(2*StretchRatio(2))), ...
        floor(ImW/2 + OutSize(2)/(2*StretchRatio(2)))+1, ...
        OutSize(2));
    Yq = linspace(floor(ImH/2 - OutSize(1)/(2*StretchRatio(1))), ...
        floor(ImH/2 + OutSize(1)/(2*StretchRatio(1)))+1, ...
        OutSize(1));
    % Create sampling meshes and do interpolation. Interpolation requires
    % input of class double and I want ImOut to be the same class as ImIn.
    [XXq, YYq] = meshgrid(Xq, Yq);
    ImOut = interp2(X, Y, double(ImIn), XXq, YYq);
    ImOut = cast(ImOut, 'like', ImIn);
end