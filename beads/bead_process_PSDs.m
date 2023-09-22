%% Process multiple sets of bead data sequentially
close all

% If you know which datasets to load, put them in fnames, e.g.
% fnames = {'water_b1_t40','water_b1_t50'}; % as many as needed.

% If you don't know which datasets to load, set these two lines to 
% fnames = {}; % and
% fIdxs = []; % then run the script.
% A numbered list of data folder names will get printed out along with all
% the fnames. If you don't want all the folders to be loaded, then choose
% the numbers you want and put them in fIdxs, e.g. fIdxs = 1:3.
% Then run again and it will print an expression for the fnames line that
% you can copy and paste.

fnames = {};
fIdxs = [];

% Experiment parameters
mPerPx = 0.062e-6;           % Camera pixel size calibration
ignoreDirs = {}; % Directories to ignore (ones without data)

% Processing parameters
timeRegularise = true;  % Regularise time vector (for data acquired with fast_acq_v9 onwards)

fitPoly = [1]; % Fit a polynomial to remove drift
fitPolyOrder = 1;       % Order of polynomial to be fitted
fpass = 0;              % High-pass frequency
fstop = [];%[92.9; 124.7; 154.7] + [-0.5 0.5];  % Bandstop frequencies - one range per row.

doFFT = false;           % Calculate FFT and maybe plot

% PSD averaging
nBlocking = 20;

% Data file parameters
forceRun = false;       % Try to take data from file and reuse as much as possible
saveData = true;        % Save data to file
dataSuff = '_psd';       % Suffix for filename when saving/loading

% Plotting parameters
loadPics = false;    % Load the TIF images
saveFigs = true;
showStack = false;   % Open the image data in ImageJ

%while true

% Get all the children directories in a struct
dirList = dir;
dirList = dirList([dirList.isdir]);
dirList = dirList(3:end);
for d = 1:length(ignoreDirs)
    dirList = dirList(~strcmp({dirList.name},ignoreDirs{d}));
end
% Check cropTs has been made with the right number of elements. Throws an
% error and gives an empty list if not.
% checkCropTs(cropTs, dirList);

% Output names
tmp = {dirList.name}'; for idx = 1:length(tmp); tmp{idx} = sprintf('%i: %s', idx, tmp{idx}); end
disp(tmp)

if isempty(fnames) 
    if isempty(fIdxs); fIdxs = 1:length(dirList); end
    fprintf('\nfnames = {')
    fprintf('''%s''',dirList(fIdxs(1)).name)
    for idx = 2:length(fIdxs)
        fprintf(',''%s''', dirList(fIdxs(idx)).name)
    end
    fprintf('};\n')
else
    psds = cell(size(fnames));
    msds = cell(size(fnames));
    stiffs = zeros(2,length(fnames));
    for fI = 1:length(fnames)
        f = fnames{fI};
        
        dataFile = [dirList(fI).name '_processed' dataSuff '.mat'];
        if forceRun || ~exist(dataFile, 'file')
            data = bead_loadData(f, false);
            
            data.opts.timeRegularisation = true;
            data.opts.bandstop = fstop;
            data.mPerPx = mPerPx;
            
            data = bead_preProcessCentres(data);
            data = bead_filter_bandstop(data);
            
            data = bead_PSD(data,'nBlocking',nBlocking);
            fh = gcf;
            ax = false(size(fh.Children));
            for idx = 1:length(fh.Children); if isa(fh.Children(idx),'matlab.graphics.axis.Axes') %#ok<ALIGN>
                    ax(idx) = true;
            end;        end
            ax = fh.Children(ax);
            title(ax(2), sprintf('f_E = %i Hz, V_{pp} = %1g',data.opts.Vfreq, data.opts.Vpp))
            
             if saveFigs
                 saveas(fh, [data.fName '_PSD.png'])
                 print(fh, '-painters','-dsvg',[data.fName '_PSD.svg'])
             end
        
            data = bead_normMSD(data);
            
            fh = gcf;
            if saveFigs
                 saveas(fh, [data.fName '_MSD.png'])
                 print(fh, '-painters','-dsvg',[data.fName '_MSD.svg'])
            end
             
        else
            tmp = load(dataFile, 'data');
            data = tmp.data;
            clear tmp
            data.opts.forceRun = forceRun;
        end
    
        
        
        psds{fI} = data.pro.psd;
        msds{fI} = data.pro.amsdObj;
        stiffs(:,fI) = 2e12 * kBT(data.opts.Temp) ./ mmsd(data.pro.amsdObj, [0.1 1]);
        
        if saveData
            save(dataFile, 'data')
        end
    end
end