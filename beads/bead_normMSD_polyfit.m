function data = bead_normMSD_polyfit(data, direction, offset, varargin)
%% data = bead_normMSD_polyfit(data, direction, offset, [num_t, doPlots])
% Take data object and use x, y or both from raw data, returning calculated
% MSD and the MSD object in complete data struct. Uses Tinevez's
% msdanalyzer for the bulk of the work

% Parse the inputs - needs to be after centres is defined (above)
doPlots = true;
if nargin == 3
    num_t = data.nPoints;
elseif nargin >= 4
    num_t = varargin{1};
    if nargin == 5
        doPlots = varargin{2};
    end
else
    error('Wrong number of input arguments')
end

% A little bit hacky - if both directions are wanted, stick them together
% and replicate offset array. This means offset is the same for both
% dimensions
nPerDir = length(offset);
if strcmp(direction, 'x') || strcmp(direction, 'y')
    centres = data.raw.([direction 'CentresPx']) * data.mPerPx;
    legCell = repmat({direction},nPerDir,1);
else
    centres = [data.raw.xCentresPx data.raw.yCentresPx] * data.mPerPx;
    timeVec = [data.raw.timeVecMs data.raw.timeVecMs];
    offset = [offset (offset + data.nPoints)];
    legCell = [repmat({'X'},1,nPerDir) repmat({'Y'},1,nPerDir)];
end

% Dimension order for polynomial fitting
dims = [1, 3, 2];

if data.opts.forceRun || (~isfield(data.pro, 'amsdObj') && ~isfield(data.pro, [direction(1) 'msdObj']))
    % Prepare data to go into msdanalyzer
    tracks = cell(length(offset),1);
    
    for idx = 1:length(offset)
        % Crop data, then fit polynomial of order set by processing script
        centresCrop = centres(offset(idx) : offset(idx) + num_t - 1);
        [~, centresCrop, ~] = func_thermal_rm(1:length(centresCrop), ...
            permute(centresCrop, dims), data.opts.pOrder, 1, length(centresCrop));
        centresCrop = ipermute(centresCrop, dims);
        
        tracks{idx} = [1e-3 .* timeVec(offset(idx) : offset(idx) + num_t - 1)' ...
            1e6 .* centresCrop'];
    end
    
    % Make an msdanalyzer and use it
    msd = msdanalyzer(1, 'um', 's');
    msd = msd.addAll(tracks);
    msd = msd.computeMSD;
    
    % Calculate normalized MSD
    MSDs = cat(3,msd.msd{:});
    dTs = squeeze(MSDs(:, 1, :));
    MSDs = squeeze(MSDs(:, 2, :));
    Xs = cat(3, msd.tracks{:});
    Xs = squeeze(Xs(:,2,:));
    Xvars = var(Xs);
    MSDnorm = 0.5*MSDs./Xvars;
    
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
    MSDnorm = data.pro.([savedDirection 'MSDnorm']);
    dTs = MSDnorm(:,1);
    MSDnorm = MSDnorm(:,2);
    tracks = data.pro.([savedDirection 'msdObj']).tracks;
end


if doPlots
    fh = figure;
    fh.Name = data.fName;
    clf
    hold on
    
    if strcmp( direction(1) , 'a')
        legCell = {[]};
        for Idx = 1:length(tracks)/2
            legCell{Idx} = [num2str(round(tracks{Idx}(1,1))) 's - ' num2str(round(tracks{Idx}(end,1))) 's'];
        end
        
        subplot(2,1,1);
        
        loglog(dTs(:,1:end/2), MSDnorm(:,1:end/2), 'LineWidth',2);
        xlim([1e-3 diff(tracks{1}([1 end/2],1))])
        xlabel('Delay (s)')
        ylabel('Normalized MSD')
        title({'Normalized mean square X displacement',['after polynomial order ' num2str(data.opts.pOrder) ' fitting'], ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
        legend(legCell,'Location','best')
        
        subplot(2,1,2);

        loglog(dTs(:,1+end/2:end), MSDnorm(:,1+end/2:end), 'LineWidth',2);
        xlim([1e-3 diff(tracks{1}([1 end/2],1))])
        xlabel('Delay (s)')
        ylabel('Normalized MSD')
        title({'Normalized mean square Y displacement',['after polynomial order ' num2str(data.opts.pOrder) ' fitting'], ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
        legend(legCell,'Location','best')
    else
        loglog(dTs, MSDnorm,'LineWidth',2)
        
        fh.Children.XAxis.Scale = 'log';
        fh.Children.YAxis.Scale = 'log';
        
        xlim([1e-3 tracks{1}(end/2,1)])
        
        xlabel('Delay (s)')
        ylabel('Normalized MSD')
        title({'Normalized mean square displacement',['after polynomial order ' num2str(data.opts.pOrder) ' fitting'], ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
        legend(legCell, 'Location','best')
    end
    % legend('20 - 40s','40 - 60s', '3m00s - 3m20s','3m20s - 3m40s','32m00s - 32m20s', '32m20s - 32m40s','Location','best')
end