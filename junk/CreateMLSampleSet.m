%% Take some images out to create training dataset
% Output images directory
SaveDir = '~/AI_Optrap/Samples/Set6/';
% How many images to take from each dataset
FramesPerSet = 100;
% Dummy class - images in this class are all 0s. 
% 'stretched' or 'relaxed'. Anything else will be ignored/treated as none
DummyClass = 'none';

% Which datasets to take images from
Cells = {'LS174T', 'HL60', 'MV411'};
DSets = {{'normoxia','hypoxia'},{'normoxia','with_drugs'},{'normoxia','with_drugs'}};
Nums = {{1:20, 1:20},{1:30, 1:20},{1:30, 1:20}};

% Calculate how many images there will be
AllNums = [Nums{:}];
NumSets = length([AllNums{:}]);
NumIms = FramesPerSet * NumSets;

SelectedIms = cell(NumIms, 2);
Count = 1;
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
                
                % Round up output size to nearest power of 2
                %OutSize = repmat(2^ceil(log2(max(size(Imstack{1}{1,1})))),1,2);
                % Output size = image size
                OutSize = size(Imstack{1}{1,1});
                % Output size = 512 * 512
                %OutSize = min([512, 512],[ImInfo.ImH, ImInfo.ImW]);
                % Output size = 256 * 256
%                 OutSize = [256, 256];
                
                SelectedIms{Count,1} = repmat(uint8(mean(Imstack{1}{ImInfo.FrNum},'all')),OutSize(1),OutSize(2));
                SelectedIms{Count,1}(1:OutSize(1), 1:OutSize(2)) = Imstack{1}{ImInfo.FrNum,1}... Need to index to take the centre of the image
                    (floor((ImInfo.ImH-OutSize(1))/2)+1:floor((ImInfo.ImH+OutSize(1))/2),... 
                    floor((ImInfo.ImW-OutSize(2))/2)+1:floor((ImInfo.ImW+OutSize(2))/2));
                SelectedIms{Count,2} = ImInfo;
                Count = Count + 1;
            end
        end
    end
end
toc(StartTime)
%% Augment
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

for Imidx = 1:Count-1
    II = SelectedIms{Imidx,2};
    if II.FrNum > 500
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
    FName = strjoin({II.CellType, II.Set, II.Num, num2str(II.FrNum)},'_');
%     if strcmp(DummyClass,ThisImClass)
        imwrite(uint8(SelectedIms{Imidx,1}),[Dir FName '.png'])
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

