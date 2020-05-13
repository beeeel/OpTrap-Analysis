%% Take some images out to create training dataset
% Output images directory
SaveDir = '~/png/ML/Samples/';
% How many images to take from each dataset
FramesPerSet = 100;

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
            LoadImstackInfoMeta(CellType, DSet, NumStr)
%             BlockSize = meta.N_Frames./FramesPerSet;
%             Offset = randi(BlockSize);            
%             for frame = 0:FramesPerSet-1
            for frame = 1:FramesPerSet
                if frame <= FramesPerSet/2
                    Offset = 0;
                elseif meta.N_Frames == 1000
                    Offset = 850 - FramesPerSet/2;
                elseif meta.N_Frames == 2000
                    Offset = 1850 - FramesPerSet/2;
                else
                    error([meta.N_Frames ' frames in this set: ' meta.filepath])
                end

                ImInfo = struct;
                ImInfo.CellType = CellType;
                ImInfo.Set = DSet;
                ImInfo.Num = NumStr;
                ImInfo.FrNum = frame + Offset; %frame * BlockSize + Offset;
                ImInfo.Radius = info(ImInfo.FrNum).radius;
                ImInfo.Centre = info(ImInfo.FrNum).centres;
                ImInfo.ImW = size(Imstack{1}{ImInfo.FrNum},2);
                ImInfo.ImH = size(Imstack{1}{ImInfo.FrNum},1);
                
                SelectedIms{Count,1} = zeros(2^ceil(log2(max(size(Imstack{1}{1,1})))));
                SelectedIms{Count,1}(1:ImInfo.ImH, 1:ImInfo.ImW) = Imstack{1}{ImInfo.FrNum,1};
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
%% Save these files as .pgm (portable gray map)
LastNum = '';
for Imidx = 1:Count-1
    II = SelectedIms{Imidx,2};
    if II.FrNum > 500
        Dir = [SaveDir 'stretched/'];
    else
        Dir = [SaveDir 'relaxed/'];
    end
    FName = strjoin({II.CellType, II.Set, II.Num, num2str(II.FrNum)},'_');
    imwrite(SelectedIms{Imidx,1},[Dir FName '.pgm'])
    [FID, msg] = fopen([Dir FName '.txt'],'w');
    fprintf(FID,'CellType\t%s\tSet\t%s\tSetNum\t%s\tFrame\t%i\tRadius\t%g\tCentre(x)\t%g\tCentre(y)\t%g',II.CellType,II.Set,II.Num,II.FrNum,II.Radius,II.Centre(1),II.Centre(2));
    fclose(FID);
    if ~strcmp(II.Num,LastNum)
        disp(['Started files for ' II.CellType ' ' II.Set ' ' II.Num])
    end
    LastNum = II.Num;
end

