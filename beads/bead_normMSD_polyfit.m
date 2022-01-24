function data = bead_normMSD_polyfit(data, direction, offset, varargin)
%% data = bead_normMSD_polyfit(data, direction, offset, [num_t, doPlots, useRaw, centresRow, doNorm, errorBars, useField])
% Take data object and use x, y or both from raw data, returning calculated
% MSD and the MSD object in complete data struct. Uses Tinevez's
% msdanalyzer for the bulk of the work

% Parse the inputs
doPlots = true;
useRaw = false;
% If working with 1bead data, use 1 row, for 2bead data, take 2.
centresRow = 1;
if strcmp(data.raw.suffixes{1}, 'l') && strcmp(data.raw.suffixes{2}, 'r')
    centresRow = [1 2];
end

if length(data.opts.cropT) == 2
    cropT = data.opts.cropT;
else
    cropT = [1 length(data.raw.timeVecMs)];
end

num_t = min(data.nPoints, diff(cropT)+1);
doNorm = true;
errorBars = false;
useField = '';
if nargin >= 4 && ~isempty(varargin{1})
    num_t = varargin{1};
end
if nargin >= 5 && ~isempty(varargin{2})
    doPlots = varargin{2};
end
if nargin >= 6 && ~isempty(varargin{3})
    useRaw = varargin{3};
end
if nargin >= 7 && ~isempty(varargin{4})
    centresRow = varargin{4};
end
if nargin >= 8 && ~isempty(varargin{5})
    doNorm = logical(varargin{5});
end
if nargin >= 9 && ~isempty(varargin{6})
    errorBars = logical(varargin{6});
end
if nargin >= 10 && ~isempty(varargin{7})
    useField = varargin{7};
end
if nargin > 10
    warning('Wrong number of input arguments')
    warning(['Ignoring ' num2str(nargin-10) ' arguments'])
end

% A little bit hacky - if both directions are wanted, stick them together
% and replicate offset array. This means offset is the same for both
% dimensions

nPerDir = length(offset);
if isfield(data.opts, 'UseField') && ~useRaw
    % If there's a field we should use, and we're not forced to use raw
    useField = get_Fn;
    
    if isfield(data.opts, 'cropTHPval')
        cropTHPval = data.opts.cropTHPval;
    else
        cropTHPval = 0;
    end
    
    cropTHP = [cropTHPval+1, diff(cropT) + 1 - cropTHPval];

    tmp = data.raw.timeVecMs(cropT(1):cropT(2));
    tmp = tmp(cropTHP(1):cropTHP(2));
    
    if any(strcmp(direction, {'x', 'y'}))
        centres = data.pro.([direction useField])(centresRow,cropT(1):cropT(2));
        timeVec = tmp;
        legCell = repmat({direction},nPerDir,1);
        filtStr = ['using field ' useField];
        % Hacky af: set num_t to min of actual num_t and previous value
        num_t = min(num_t, size(data.pro.(['x' useField]),2));
    else
        if size(data.pro.(['x' useField{1}]),2) > cropT(2)
            centres = [data.pro.(['x' useField{1}])(centresRow,cropT(1):cropT(2)) data.pro.(['y' useField{2}])(centresRow,cropT(1):cropT(2))];
        else
            centres = [data.pro.(['x' useField{1}])(centresRow,:) data.pro.(['y' useField{2}])(centresRow,:)];
            warning(['Didn''t apply cropT because it looks like data.pro.x' useField{1} ' has already been cropped']) 
        end
        timeVec = [tmp tmp];
        offset = [offset (offset + num_t)];
        legCell = [repmat({'X'},1,nPerDir) repmat({'Y'},1,nPerDir)];
        % Hacky af: set num_t to min of actual num_t and previous value
        num_t = min(num_t, size(data.pro.(['x' useField{1}]),2));
        if strcmp(useField{1}, useField{2})
            filtStr = ['using field ' useField{1}];
        else
            filtStr = sprintf('using fields %s and %s',useField{1}, useField{2});
        end
    end
else
    % Otherwise, if told to, use raw
    if (strcmp(direction, 'x') || strcmp(direction, 'y'))
        centres = data.raw.([direction 'CentresPx'])(centresRow,cropT(1):cropT(2)) * data.mPerPx;
        timeVec = data.raw.timeVecMs(cropT(1):cropT(2));
        legCell = repmat({direction},nPerDir,1);
    else
        centres = [data.raw.xCentresPx(centresRow,cropT(1):cropT(2)) ...
            data.raw.yCentresPx(centresRow,cropT(1):cropT(2))] * data.mPerPx;
        timeVec = [data.raw.timeVecMs(cropT(1):cropT(2)) ...
            data.raw.timeVecMs(cropT(1):cropT(2))];
        
        offset = [offset (offset + num_t)];
        legCell = [repmat({'X'},1,nPerDir) repmat({'Y'},1,nPerDir)];
    end
    filtStr = ['unfiltered'];
end
clear tmp


if data.opts.forceRun || (~isfield(data.pro, 'amsdObj') && ~isfield(data.pro, [direction(1) 'msdObj']))
    
    % % Dimension order for polynomial fitting
    % dims = [1, 3, 2];
    
    % Store which row of centres we're using
    if ~isfield(data.raw,'suffixes')
        warning('Data does not contain suffixes field, cannot determine centroid method')
    elseif isfield(data.opts, 'xHPSuffix')
        data.opts.msdSuffix = data.opts.xHPSuffix;
    elseif isfield(data.raw, 'suffixes')
        data.opts.msdSuffix = data.raw.suffixes(centresRow);
    else
        disp('data.opts does not contain xHPSuffix, and centres does not contain enough rows to use specified centres row')
        error('Cannot determine which centroid method has been used')
    end
    
    % Prepare data to go into msdanalyzer
    % When taking multiple centresRows, i.e.: 2 beads, we want all the data
    % for one bead, followed by all the data for the second.
    tracks = cell(length(offset), length(centresRow));
    for idx = 1:length(offset)
        % Crop data, then demean. (Don't fit poly because either this has
        % already been done, or it's not wanted)
        centresCrop = centres(:, offset(idx) : offset(idx) + num_t - 1);
        centresCrop = centresCrop - mean(centresCrop,2);
%         
%         [~, centresCrop, ~] = func_thermal_rm(1:length(centresCrop), ...
%             permute(centresCrop, dims), data.opts.pOrder*(~useRaw), 1, length(centresCrop));
%         centresCrop = ipermute(centresCrop, dims);

        for row = 1:length(centresRow)
            tracks{idx, row} = [1e-3 .* timeVec(offset(idx) : offset(idx) + num_t - 1)' ...
                1e6 .* centresCrop(row,:)'];
        end
    end
    
    tracks = reshape(tracks, [], 1);
    
    % Make an msdanalyzer and use it
    msd = msdanalyzer(1, 'um', 's','log');
    msd = msd.addAll(tracks);
    msd = msd.computeMSD;
    
    % Calculate normalized MSD
    MSDs = cat(3,msd.msd{:});
    dTs = squeeze(MSDs(:, 1, :));
    if errorBars
        msdStd = squeeze(MSDs(:,3,:));
    else
        msdStd = [];
    end
    MSDs = squeeze(MSDs(:, 2, :));
    Xs = cat(3, msd.tracks{:});
    Xs = squeeze(Xs(:,2,:));
    % Normalize unless told not to
    if doNorm
        Xvars = var(Xs);
        MSDnorm = 0.5*MSDs./Xvars;
    else
        MSDnorm = MSDs;
    end    
    
    % Put the nomalized MSD and msdanalyzer object into data struct
    data.pro.([direction(1) 'MSDnorm']) = [dTs MSDnorm];
    data.pro.([direction(1) 'msdObj']) = msd;
elseif doPlots
    % Figure out what the saved data is called
    savedDirection = direction(1);
    if isfield(data.pro, 'amsdObj')
        savedDirection = 'a';
    elseif ~isfield(data.pro, [direction(1) 'msdObj'])
        error('Uuugghhh you said this wouldn''t happen but it has')
    end
    
    % Extract the saved data to be used for plots below
    msd = data.pro.([savedDirection 'msdObj']);
    MSDs = cat(3,msd.msd{:});
    dTs = squeeze(MSDs(:, 1, :));
    if errorBars
        msdStd = squeeze(MSDs(:,3,:));
    else
        msdStd = [];
    end
    MSDs = squeeze(MSDs(:, 2, :));
    Xs = cat(3, msd.tracks{:});
    Xs = squeeze(Xs(:,2,:));
    % Normalize unless told not to
    if doNorm
        Xvars = var(Xs);
        MSDnorm = 0.5*MSDs./Xvars;
    else
        MSDnorm = MSDs;
    end    
    tracks = msd.tracks;
    
end


if doPlots
    do_plots;
end

%%
    function fn = get_Fn
        
        if any(strcmp(direction, {'x', 'y'}))
            if isfield(data.pro, [direction useField])
                fn = useField;
            elseif isfield(data.pro, [direction 'CentresM'])
                warning('Falling back to using CentresM because useField failed to resolve')
                fn = 'CentresM';
            else
                error('Could not decide a field to use, tried %s and CentresM', useField)
            end
        else
            fn = cell(2,1);
            xy = 'xy';
            for i = 1:2
                d = xy(i);
                if isfield(data.pro, [d useField])
                    fn{i} = useField;
                elseif isfield(data.pro, [d 'CentresM'])
                    fn{i} = 'CentresM';
                else
                    error('Could not decide a field to use, tried %s and CentresM', useField)
                end
            end
        end
        
    end

    function do_plots
        fh = figure;
        fh.Name = data.fName;
        clf
        hold on
        
        if strcmp( direction(1) , 'a')
            legCell = {[]};
            for Idx = 1:length(tracks)/2
                legCell{Idx} = [num2str(round(tracks{Idx}(1,1))) 's - ' num2str(round(tracks{Idx}(end,1))) 's'];
            end
            
            ax = subplot(2,1,1);
            if ~isempty(msdStd)
                errorbar(dTs(:,1:end/2), MSDnorm(:,1:end/2), msdStd(:,1), 'LineWidth',2);
            else
                plot(dTs(:,1:end/2), MSDnorm(:,1:end/2), 'LineWidth',2);
            end
            ax.XScale = 'log';
            ax.YScale = 'log';
            % Set X lim to show shortest delay up to half total time
            xlim([diff(tracks{1}([1 2])) diff(tracks{1}([1 round(size(tracks{1},1)/2)],1))])
            xlabel('Delay (s)')
            if doNorm
                ylabel('Normalized MSD')
                title({'Normalized mean square X displacement', filtStr , ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            else
                ylabel('MSD (μm^2)')
                title({'Mean square X displacement', filtStr , ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            end
            legend(legCell,'Location','best')
            
            ax = subplot(2,1,2);
            
            if ~isempty(msdStd)
                errorbar(dTs(:,1+end/2:end), MSDnorm(:,1+end/2:end), msdStd(:,2), 'LineWidth',2);
            else
                plot(dTs(:,1+end/2:end), MSDnorm(:,1+end/2:end),  'LineWidth',2);
            end
            ax.XScale = 'log';
            ax.YScale = 'log';
            % Set X lim to show shortest delay up to half total time
            xlim([diff(tracks{1}([1 2])) diff(tracks{1}([1 round(size(tracks{1},1)/2)],1))])
            xlabel('Delay (s)')
            if doNorm
                ylabel('Normalized MSD')
                title({'Normalized mean square Y displacement',filtStr, ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            else
                ylabel('MSD (μm^2)')
                title({'Mean square Y displacement', filtStr , ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            end
            legend(legCell,'Location','best')
        else
            ax = gca;
            errorbar(dTs, MSDnorm, msdStd, 'LineWidth',2)
            ax.XScale = 'log';
            ax.YScale = 'log';
            
            fh.Children.XAxis.Scale = 'log';
            fh.Children.YAxis.Scale = 'log';
            
            % Set X lim to show shortest delay up to half total time
            xlim([diff(tracks{1}([1 2])) diff(tracks{1}([1 round(size(tracks{1},1)/2)],1))])
            
            xlabel('Delay (s)')
            if doNorm
                ylabel('Normalized MSD')
                title({'Normalized mean square displacement',filtStr, ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            else
                ylabel('MSD (μm^2)')
                title({'Mean square displacement',filtStr, ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            end
            legend(legCell, 'Location','best')
        end
        % legend('20 - 40s','40 - 60s', '3m00s - 3m20s','3m20s - 3m40s','32m00s - 32m20s', '32m20s - 32m40s','Location','best')
    end

end

% % Some crap that I don't use anymore
%{ 
%   elseif isfield(data.opts,'fpass')
%     % Else use highpass if available
% 
%     % Hacky af: set num_t to min of actual num_t and previous value
%     num_t = min(num_t, size(data.pro.xCentresHP,2));
%     
%     cropTHPval = data.opts.cropTHPval;
%     cropTHP = [cropTHPval+1, diff(cropT) + 1 - cropTHPval];
% 
%     tmp = data.raw.timeVecMs(cropT(1):cropT(2));
%     tmp = tmp(cropTHP(1):cropTHP(2));
%     if (strcmp(direction, 'x') || strcmp(direction, 'y'))
%         centres = data.pro.([direction 'CentresHP']);
%         timeVec = tmp;
%         legCell = repmat({direction},nPerDir,1);
%     else
%         centres = [data.pro.xCentresHP data.pro.yCentresHP];
%         timeVec = [tmp tmp];
%         offset = [offset (offset + size(data.pro.xCentresHP,2))];
%         legCell = [repmat({'X'},1,nPerDir) repmat({'Y'},1,nPerDir)];
%     end
%     filtStr = ['after ' num2str(data.opts.fpass) 'Hz highpass filtering '];
% elseif isfield(data.pro,[direction 'CentresM']) || (strcmp(direction(1), 'a') && min(isfield(data.pro,{'xCentresM','yCentresM'})))
%     % If data has not been highpass filtered, use polyfiltered
%     if (strcmp(direction, 'x') || strcmp(direction, 'y'))
%         centres = data.pro.([direction 'CentresM'])(centresRow,cropT(1):cropT(2));
%         timeVec = data.raw.timeVecMs(cropT(1):cropT(2));
%         legCell = repmat({direction},nPerDir,1);
%     else
%         centres = [data.pro.xCentresM(centresRow,cropT(1):cropT(2)) ...
%             data.pro.yCentresM(centresRow,cropT(1):cropT(2))];
%         timeVec = [data.raw.timeVecMs(cropT(1):cropT(2)) ...
%             data.raw.timeVecMs(cropT(1):cropT(2))];
%         
%         nP = min(data.nPoints, diff(cropT)+1);
%         offset = [offset (offset + nP)];
%         legCell = [repmat({'X'},1,nPerDir) repmat({'Y'},1,nPerDir)];
%     end
%     filtStr = ['after polynomial order ' num2str(data.opts.pOrder) ' fitting'];
% end
%}