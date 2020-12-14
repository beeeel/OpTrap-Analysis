%% First attempt at bead analysis
% Loads images from B09, calculates centre of mass and stiffness for each
% set. Attempt to calculate accurate frame timing below.
%%
% Load the metadata, extract frame time, then load image data

% All the bead files I collected in the first day
AllBeads = {...'2020_10_19/bead_1m-ims_30pc-laser_2/Untitled_1/Untitled_1_MMStack_Default',...
    '2020_10_19/bead_1m-ims_30pc-laser_1/Untitled_1/Untitled_1_MMStack_Default',...
    '2020_10_19/bead_10k-ims_30pc-laser_1/bead_10k-ims_30pc-laser_1_MMStack_Default',...
    '2020_10_19/bead_25k-ims_30pc-laser_1/bead_25k-ims_30pc-laser_1_MMStack_Default',...
    '2020_10_19/bead_28k-ims_30pc-laser_1/bead_28k-ims_30pc-laser_1_MMStack_Default'};

DataDir = '~/Documents/data/OpTrap/bead_analysis/';
Stiffness = zeros(2, length(AllBeads));
AllFrameTimes = zeros(1, length(AllBeads));

CofM = cell(1,length(AllBeads));
%% Load the data from tif and metadata from txt and extract stuff + process
Set = 1:length(AllBeads);

for idx = Set
    FileBase = AllBeads{idx};
    tic
    metadata = fileread([DataDir FileBase '_metadata.txt']);
    metadata = jsondecode(metadata);
    
    %{
    % Frame time in ms
    FrameTime = str2double(metadata.FrameKey_0_0_0.Camera_1_ActualInterval_ms);
    AllFrameTimes(idx) = FrameTime;
    clear metadata Imstack
    disp(['Finished parsing metadata at ' num2str(toc) 's. Starting Imstack'])
    Imstack = bfopen([DataDir FileBase '.ome.tif']);
    disp(['Finished reading Imstack at ' num2str(toc) 's. Starting analysis'])
    %%
    Ims = double(cat(3,Imstack{1}{:,1}));
    FileName = strsplit(AllBeads{idx},'/');
    save(['~/Documents/data/OpTrap/bead_analysis/MatFiles/' FileName{1} '_' FileName{2}],'Ims','FrameTime','-v7.3');
    
    Sz = size(Ims);
    Ims = reshape(Ims,[],Sz(3));
    Ims = normalize(Ims,'center');
    Ims = reshape(Ims.^2,Sz(1), Sz(2), Sz(3));
    [X, Y] = meshgrid(1:size(Ims,2),1:size(Ims,1));
    % Calculate centre of mass - sum weighted by pixel location, normalized to
    % image total.
    CofM{idx} = squeeze([sum(Ims.*X,[1 2]) sum(Ims.*Y,[1 2])]./sum(Ims,[1 2]));
    clear Ims
    %% Stiffness = Kb * T * var(x)
    Kb = 1.38064852e-23; % Boltzmann constant
    T = 273 + 20; % Assume trap at room temperature - maybe I should model this?
    % From Sarshar et al 2014 eq 4. Need to convert CofM from px to m, so this
    % gives stiffness in N/m.
    Stiffness(:,idx) = Kb .* T ./ var(CofM{idx}*0.07e-6,0,2);
    %}
end
%% Extract frame timing and intervals from the metadata
FNames = fieldnames(metadata);
FTimes = zeros(1,length(FNames)-1);
FInts = zeros(1,length(FNames)-1);
for idx = 2:length(FNames)
    FTimes(idx-1) = datenum(metadata.(FNames{idx}).ReceivedTime,'yyyy-mm-dd HH:MM:SS.FFF');
    FInts(idx-1) = str2double(metadata.(FNames{idx}).Camera_1_ActualInterval_ms);
end
%% Show some stuff
% Predicted time is "ActualIntervalMs" * number of intervals elapsed
TPred = FInts(1) .*(0:length(FNames)-2);
% Elapsed time is received time (this frame) - received time (first frame)
TElapsedDays = (FTimes-FTimes(1));
% Convert to ms
TElapsedMs = TElapsedDays * 24 * 60 * 60 * 1000;
% Difference is predicted - elapsed
TDiff = TPred - TElapsedMs;
% Actual (mean) interval is total elapsed / number of frames
IntActualMs = 24 * 60 * 60 * 1e3 * (FTimes(end) - FTimes(1)) ./ (length(FNames)-1);

figure(3)
subplot(2,1,1)
plot(1:length(FNames)-1,FInts)
subplot(2,1,2)
%% Show an image with the centre of mass on it
Frame = 10311;
figure( 1)
imagesc(Imstack{1}{Frame, 1})
hold on
axis image
plot(CofM{end}(1,Frame), CofM{end}(2,Frame),'kx')
%% Scatter plot and position-time plots for each set
for SetN = 1:length(CofM)
    figure(SetN)
    clf
    plot(CofM{SetN}(1,:).*0.07,0.07.*CofM{SetN}(2,:),'.')
    hold on
    plot(0.07.*mean(CofM{SetN}(1,:)), 0.07.*mean(CofM{SetN}(2,:)),'rx')
    title(['Scatter plot of position from ' num2str(length(CofM{SetN})) ' frames'])
    xlabel('X (\mu m)')
    ylabel('Y (\mu m)')
    legend('Measured position','Mean position')
    
    FileName = strsplit(AllBeads{SetN},'/');
%     saveas(gcf,['~/Documents/data/OpTrap/bead_analysis/MatFiles/' FileName{1} '_' FileName{2} '.png'])
    
    figure(SetN + length(CofM))
    clf
    subplot(2,1,1)
    plot(0.07.*CofM{SetN}(1,:))
    title(['X position from ' num2str(length(CofM{SetN})) ' frames'])
    xlabel('Frame index')
    ylabel('Position (\mu m)')
    subplot(2,1,2)
    plot(0.07.*CofM{SetN}(2,:))
    title(['Y position from ' num2str(length(CofM{SetN})) ' frames'])
    xlabel('Frame index')
    ylabel('Position (\mu m)')
end
%% Load the .mat files to look at the videos
FileNum = 2;
FileName = strsplit(AllBeads{FileNum},'/');
data = load(['~/Documents/data/OpTrap/bead_analysis/MatFiles/' FileName{1} '_' FileName{2}]); 

Frs = 1:1000:size(data.Ims,3);
figure
for frame = Frs
    imagesc(data.Ims(:,:,frame))
    pause(0.25)
end
%% Write hashtable to text file
metaFileName = ['~/Documents/data/Zeiss/2020-11-10-2um-beads-broken-XT-files/' fileName(1:end-3) 'txt'];
system(['touch ' metaFileName]);
hashTable = Metastack{2};
hashTableChar = char(hashTable.toString());
% dlmwrite(metaFileName,hashTableChar,'delimiter','');
keysCell = strsplit(hashTableChar(2:end-1),',')';
for idx = 1:size(keysCell,1)
    split = strsplit(keysCell{idx},'=');
    keysCell{idx,1} = split{1};
    keysCell{idx,2} = split{2};
end
 