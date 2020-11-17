%% First attempt at bead analysis
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
%% Load some confocal bead data
fileName = 'laser-25-100kfr_few_px.czi';
start = tic;
Imstack = bfopen(['~/Documents/data/Zeiss/2020-11-10-2um-beads-repaired-XY-files/' fileName]);
imTime = toc(start);
Metastack = bfopen(['~/Documents/data/Zeiss/2020-11-10-2um-beads-broken-XT-files/' fileName]);
metaTime = toc(start) - imTime;
%%
% Metadata handling
omeMeta = Metastack{4};
hashTable = Metastack{2};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeT = omeMeta.getPixelsSizeT(0).getValue(); % stack height, lines

voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.METER); % in µm
voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double
voxelSizeTdouble = str2double(hashTable.get('Global Information|Image|Channel|LaserScanInfo|LineTime #1'));

% Image processing
if size(Imstack{1},1) == 1
    Ims = Imstack{1}{1,1};
else
    warning('untested line of code ahead. Check Ims variable is what it should be');
    Ims = cat(2,Imstack{1}{:,1});
end

% Crop to just the bead
cropX = [50, 70];
cropT = [1 8.5e4];
ImsC = zeros(cropT(2)/2, cropX(2)-cropX(1)+1);
ImsC(:,:,1) = Ims(cropT(1):2:cropT(2), cropX(1):cropX(2));
ImsC(:,:,2) = Ims(cropT(1)+1:2:cropT(2), cropX(1):cropX(2));
% Calculate centres in units of m (determined by unit set getting
% voxelSizeX)
ImsT = double(permute(ImsC,[2 1 3]));
X = voxelSizeXdouble.*(1:size(ImsT,1))';
Centres = sum(X.*ImsT)./sum(ImsT);
centresIL(1:2:diff(cropT)+1) = Centres(:,:,1);
centresIL(2:2:diff(cropT)+1) = Centres(:,:,2);

% Calculate bead widths in units of um
thresholded = abs(ImsT) > 50;
% If D(i+1) > D(i) you're in a new block of changing position
shiftedD = zeros(size(thresholded),'logical');
shiftedD(2:end,:) = thresholded(1:end-1,:);
beadStart = shiftedD < thresholded;
% If D(i-1) > D(i) you're in a new block of equilibrium
shiftedU = zeros(size(thresholded),'logical');
shiftedU(1:end-1,:) = thresholded(2:end,:);
beadEnd = shiftedU < thresholded;
% max returns the first true from each column
[~, beadsStartX]= max(beadStart);
[~, beadsEndX] = max(flipud(beadEnd));
% Blocks ready for output
beadsX = [beadsStartX; beadsEndX];
beadWidthsX = diff(beadsX);

% Calculate stiffness in units N/m (or uN/um)
Kb = 1.38064852e-23; % Boltzmann constant
T = 273 + 20; % Assume trap at room temperature - maybe I should model this?
% From Sarshar et al 2014 eq 4. Need to convert CofM from px to m, so this
% gives stiffness in N/m.
varCentres = var(Centres,0,2);
varStart = var(beadsStartX,0,2);
varEnd = var(beadsEndX,0,2);
Stiffness = Kb .* T ./ sqrt(0.5 * (varCentres(1).^2 + varCentres(2).^2));
StiffnessL = Kb .* T ./ sqrt(0.5 * (varStart(1).^2 + varStart(2).^2));
StiffnessR = Kb .* T ./ sqrt(0.5 * (varEnd(1).^2 + varEnd(2).^2));

% Plot stuff
figure(1)
clf
m = 3;

subplot(m,1,1) % Cropped bead images
imagesc([0 voxelSizeTdouble*diff(cropT)],1e6.*[X(1) X(end)],Ims(cropT(1):cropT(2), cropX(1):cropX(2))')
set(gca,'YDir','normal')
xlabel('Time(s)')
ylabel('Scan Axis (\mu m)')
title('Cropped data (85,000 frames)')

subplot(m,1,2) % Histogram of centres from two scan directions
histogram(reshape(Centres(:,:,1),1,[],1))
hold on
histogram(reshape(Centres(:,:,2),1,[],1))
ylabel('Count')
xlabel('Position (\mu m)')
title('Bead position distribution (bidirectional scan)')
legend('Scan direction 1','Scan direction 2')

subplot(m,1,3) % Centre position against time
plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*centresIL,'k.')
xlabel('Time(s)')
ylabel('Position (\mu m)')
xlim(cropT*voxelSizeTdouble)
title(['Bead positions (Trap stiffness ' num2str(mean(Stiffness)*1e6) 'pN/um)'])
%%
m = 2;
figure(2)
subplot(m,1,1)
clf
hold on
plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*voxelSizeXdouble.* reshape(beadsStartX,1,[],1))
plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*voxelSizeXdouble.* reshape(beadsEndX,1,[],1))
plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*reshape(Centres,1,[],1),'k.')
xlabel('Time(s)')
ylabel('Position (\mu m)')
xlim(cropT*voxelSizeTdouble)
title(['Bead positions (Trap stiffness ' num2str(mean(Stiffness)*1e6) 'pN/um)'])
legend('LHS','RHS','Centre','Location','southeast')
subplot(m,1,2)
histogram(voxelSizeXdouble*1e6*beadWidthsX)
xlim([0, 5])
title(['Measured width of 2 \mu m bead (pixel size ' num2str(voxelSizeXdouble*1e6) '\mu m)'])

%% 
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
 