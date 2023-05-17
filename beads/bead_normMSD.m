function data = bead_normMSD(data, varargin)
%% data = bead_normMSD(data, ...)
% Take data struct and calculate MSD for processed or raw data, returning
% calculated MSD and the MSD object in complete data struct. Uses Tinevez's
% msdanalyzer for the bulk of the work.
%
% Additional parameters in the form of name-value pairs. Possible options:
%   direction       - Direction used for MSD calculation. 'x','y', or 'a'
%   offset          - Offset from start. Default 1
%   numT           - Number of time points to use. Default all
%   doPlots         - Do the plots if true, else just the calculations
%   useRaw          - Use raw data instead of processed (you probably don't want this)
%   centresRow      - Which row of the centres matrix to use. Default all.
%   doNorm          - Normalize (divide) the MSD by the variance of the position
%   errorBars       - Plot errorbars (+/- 1 s.d.) on the MSD
%   useField        - Specify which processed data field to use.


% Parse the inputs
p = inputParser;

p.addRequired('data',@(x) isa(x,'struct') && isscalar(x) );

p.addParameter('forceRun',false, @(x)islogical(x))
p.addParameter('direction','a',@(x)any(strcmp(x,{'a','all','x','y'})))
p.addParameter('offset',1,@(x)validateattributes(x,{'numeric'},{'<',data.nPoints}))
p.addParameter('numT',[],@(x)validateattributes(x,{'numeric'},{'<=',data.nPoints}))
p.addParameter('doPlots',true, @(x)islogical(x))
p.addParameter('useRaw',false,@(x)islogical(x))
p.addParameter('centresRow',[],@(x)validateattributes(x,{'numeric'},{'positive','integer','<=',numel(data.raw.suffixes)}))
p.addParameter('doNorm',false,@(x)islogical(x))
p.addParameter('errorBars', false, @(x)islogical(x))
p.addParameter('useField', [], @(x)any([isfield(data.pro,x),isfield(data.raw,x)]))
p.addParameter('compliance',false, @(x)islogical(x) && isfield(data.opts,'beadDiam'))

% p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3))
% p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
% p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(data, varargin{:});

direction = p.Results.direction;
offset = p.Results.offset;
doPlots = p.Results.doPlots;
useRaw = p.Results.useRaw;
doNorm = p.Results.doNorm;
errorBars = p.Results.errorBars;
useField = p.Results.useField;
forceRun = p.Results.forceRun || data.opts.forceRun;
compliance = p.Results.compliance;

% If working with 1bead data, use 1 row, for 2bead data, take 2.
centresRow = p.Results.centresRow;
if isempty(centresRow)
    if isfield(data.opts, 'centresRow')
        centresRow = data.opts.centresRow;
    else
        centresRow = 1:length(data.raw.suffixes);
    end
end
if isfield(data.raw,'suffixes') && strcmp(data.raw.suffixes{1}, 'l') && strcmp(data.raw.suffixes{2}, 'r') && isscalar(centresRow)
    warning('Looks like you have 2 bead data, but I''m only using 1 of them')
end

cropT = data.opts.cropT;

num_t = min([p.Results.numT, data.nPoints, diff(cropT)+1]);
if ~isempty(p.Results.numT) && num_t ~= p.Results.numT
    warning('Number of time points overriden')
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

    tmp = data.pro.timeVecMs(cropT(1):cropT(2));
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
            %warning(['Didn''t apply cropT because it looks like data.pro.x' useField{1} ' has already been cropped']) 
        end
        timeVec = [tmp tmp];
        offset = [offset (offset + length(tmp))];
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
        
        offset = [offset (offset + length(tmp))];
        legCell = [repmat({'X'},1,nPerDir) repmat({'Y'},1,nPerDir)];
    end
    filtStr = ['unfiltered'];
end
clear tmp


if forceRun || ~isfield(data.pro, 'amsdObj') || ~isfield(data.pro, [direction(1) 'msdObj'])
    
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
    tracks = cell(length(offset), length(num_t), length(centresRow));
    for idx = 1:length(offset)
        for jdx = 1:length(num_t)
            % Crop data, then demean.
            centresCrop = centres(:, offset(idx) : offset(idx) + num_t(jdx) - 1);
            centresCrop = centresCrop - mean(centresCrop,2);
            
            for row = 1:length(centresRow)
                tracks{idx, jdx, row} = [1e-3 .* timeVec(offset(idx) : offset(idx) + num_t(jdx) - 1)' ...
                    1e6 .* centresCrop(row,:)'];
            end
        end
    end
    
    tracks = reshape(tracks, [], 1);
    
    % Make an msdanalyzer and use it
    msd = msdanalyzer(1, 'um', 's','log');
    msd = msd.addAll(tracks);

    if ~compliance
        msd = msd.computeMSD;
    else
        % Assume bead diameter is in microns, apply conversion for MSD from um^2 to m^2
        beadRadius = 0.5 * data.opts.beadDiam * 1e-6 * 1e-12; 
        msd = msd.computeMSD([], beadRadius);
    end
    
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
    
else
    warning('Didn''t really do anything. Does data already have an msdObj, and no forceRun set?')
end

if doPlots
    do_plots;
end

%%
    function fn = get_Fn
        
        if isempty(useField) && isfield(data.opts, 'UseField')
            useField = data.opts.UseField;
            fprintf('data.opts says use %s for MSD calculation.\n',useField)
        elseif any(strcmp(useField(1), {'x','y'}))
            useField = useField(2:end);
        end
        
        if any(strcmp(direction, {'x', 'y'}))
            if isfield(data.pro, [direction useField])
                fn = useField;
                %fprintf('Going to calculate MSD using field: %s%s\n',direction,useField)
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
                    %fprintf('Going to calculate MSD using field: %s%s\n',d,fn{i})
                elseif isfield(data.pro, [d 'CentresM'])
                    fn{i} = 'CentresM';
                    fprintf('Falling back to calculating MSD using field: %s%s\n',d,fn{i})
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
                legCell{Idx,1} = [num2str(round(tracks{Idx}(1,1))) 's - ' num2str(round(tracks{Idx}(end,1))) 's'];
            end
            
            ax = gca;%subplot(2,1,1);
            hold on
            if ~isempty(msdStd)
                h(1) = errorbar(dTs(:,1:end/2), MSDnorm(:,1:end/2), msdStd(:,1), 'LineWidth',2);
                h(2) = errorbar(dTs(:,1+end/2:end), MSDnorm(:,1+end/2:end), msdStd(:,2), 'LineWidth',2);
            else
                h = plot(dTs(:,1:end/2), MSDnorm(:,1:end/2), ...
                    dTs(:,1+end/2:end), MSDnorm(:,1+end/2:end), 'LineWidth',2);
            end

            ax.XScale = 'log';
            ax.YScale = 'log';
            % Set X lim to show shortest delay up to half total time
            xlim([diff(tracks{1}([1 2])) diff(tracks{1}([1 round(size(tracks{1},1)/2)],1))])
            xlabel('Delay (s)')
            if doNorm
                ylabel('Normalized MSD')
                title({'Normalized mean square displacement',filtStr, ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            else
                ylabel('MSD (μm^2)')
                title({'Mean square displacement', filtStr , ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
            end
            
%             x = {dTs(:,1:end/2), dTs(:,1+end/2:end)};
%             y = {MSDnorm(:,1:end/2), MSDnorm(:,1+end/2:end)};
            hold(ax, 'on')
            cols = {ax.Children.Color}; %{[0 0.447 0.741] [0.85 0.325 0.098]};
            cols = cols(end:-1:1);
            for Idx = 1:size(dTs,2)
                [dydx, tout] = msd_gradientor(dTs(:,Idx), MSDnorm(:,Idx), 'lsq', 10);
                yyaxis(ax, 'right')
                semilogx(ax, tout, dydx, '--','LineWidth',2, 'Color',cols{Idx})
            end
            plot(ax, xlim(ax), [1 1], ':','Color',0.75*[1 1 1], 'LineWidth', 3)
            ylim(ax, [0 2])
            
            %             legend(legCell,'Location','best')
            legend(h, 'X','Y')
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
        drawnow
    end

end
