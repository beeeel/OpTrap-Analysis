function data = bead_ACF(data, varargin)
%% data = bead_ACF(data, ...)
% Take data struct and calculate ACF for processed or raw data, returning
% calculated ACF in complete data struct. 
%
% Additional parameters in the form of name-value pairs. Possible options:
%   direction       - Direction used for ACF calculation. 'x','y', or 'a'
%   doPlots         - Do the plots if true, else just the calculations
%   useRaw          - Use raw data instead of processed (you probably don't want this)
%   centresRow      - Which row of the centres matrix to use. Default all.
%   doNorm          - Normalize (divide) the ACF by the variance of the position
%   useField        - Specify which processed data field to use.
%   forceRun        - Force analysis to run, even if calculation was already done for this data
%   nAvgs           - Break track into shorter section and average ACF for each section - improves SNR, default 10

% Parse the inputs
p = inputParser;

p.addRequired('data',@(x) isa(x,'struct') && isscalar(x) );

p.addParameter('forceRun',false, @(x)islogical(x))
p.addParameter('direction','a',@(x)any(strcmp(x,{'a','all','x','y'})))
p.addParameter('doPlots',true, @(x)islogical(x))
p.addParameter('doFits',isfield(data.opts,'Vfreq'), @(x)islogical(x))
p.addParameter('useRaw',false,@(x)islogical(x))
p.addParameter('centresRow',[],@(x)validateattributes(x,{'numeric'},{'positive','integer','<=',numel(data.raw.suffixes)}))
p.addParameter('doNorm',false,@(x)islogical(x))
p.addParameter('useField', [], @(x)any([isfield(data.pro,x),isfield(data.raw,x)]))
p.addParameter('nAvgs',10,@(x)validateattributes(x,{'numeric'},{'positive','integer','<=',data.nPoints}))
p.addParameter('fitCycles',20,@(x)validateattributes(x,{'numeric'},{'positive','integer'}))

% p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3))
% p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
% p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(data, varargin{:});

direction = p.Results.direction;
doPlots = p.Results.doPlots;
doFits = p.Results.doFits;
useRaw = p.Results.useRaw;
doNorm = p.Results.doNorm;
useField = p.Results.useField;
forceRun = p.Results.forceRun || data.opts.forceRun;
nAvgs = p.Results.nAvgs;
centresRow = p.Results.centresRow;
cropT = data.opts.cropT;
fitCycles = p.Results.fitCycles;

if isfield(data.pro, 'acf') && ~forceRun
    acfs = data.pro.acf(:,2:end);
    lags = data.pro.acf(:,1);
    if any(strcmp(direction, {'x', 'y'}))
        legCell = {direction};
    else
        legCell = {'X','Y'};
    end
else
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

    if isfield(data.opts, 'UseField') && ~useRaw
        % If there's a field we should use, and we're not forced to use raw
        useField = get_Fn;

        timeVec = 1e-3 * data.pro.timeVecMs(cropT(1):cropT(2));

        if any(strcmp(direction, {'x', 'y'}))
            centres = data.pro.([direction useField])(centresRow,:);
            legCell = {direction};
        else
            centres = [data.pro.(['x' useField{1}])(centresRow,:); data.pro.(['y' useField{2}])(centresRow,:)];
            legCell = {'X' 'Y'};

        end
    else
        timeVec = 1e-3 * data.raw.timeVecMs(cropT(1):cropT(2));
        % Otherwise, if told to, use raw
        if (strcmp(direction, 'x') || strcmp(direction, 'y'))
            centres = data.raw.([direction 'CentresPx'])(centresRow,cropT(1):cropT(2)) * data.mPerPx;
            legCell = {direction};
        else
            centres = [data.raw.xCentresPx(centresRow,cropT(1):cropT(2)); ...
                data.raw.yCentresPx(centresRow,cropT(1):cropT(2))] * data.mPerPx;

            legCell = {'X' 'Y'};
        end
    end

    ntau = size(centres,2)/nAvgs;
    acfs = zeros(2*ntau-1, size(centres,1), nAvgs);
    for idx = 1:size(centres,1)
        for jdx = 1:nAvgs
            jdxs = 1+(jdx-1)*ntau:jdx*ntau;
            [acfs(:,idx,jdx), lags] = xcorr(centres(idx, jdxs));
        end
    end
    acfs = mean(acfs,3);

    if doNorm
        normF = var(centres,0,2);
        acfs = acfs ./ normF;
    end

    lags = lags * diff(timeVec([1 2]));

    data.pro.acf = [reshape(lags,[],1) acfs];

    if doFits
        if isfield(data.opts,'Vfreq')
            % Here's the ACF function that you wanna fit to.
            % p = [gamma, tauc]
            fnc = @(gamma, tauc, x) 1 ./ (1+gamma.^2) .* exp(-x / tauc) + gamma.^2./(1+gamma.^2) .* cos(2*pi*data.opts.Vfreq.*x);
            ft = fittype( fnc);
            fopt = fitoptions('method','Nonlin','StartPoint',[0, 0.1],'Lower',[0 0],'Upper',[10, 10]);

            inds = find(lags >= 0 & lags <= fitCycles./data.opts.Vfreq);
        else
            fnc = @(tauc, x) exp(-x / tauc);
            ft = fittype(fnc);
            fopt = fitoptions('method','Nonlin','StartPoint',[0.1],'Lower',[0],'Upper',[10]);

            inds = find(lags >= 0 & lags <= 1);
        end

        nacf = acfs(inds,:) ./ acfs(inds(1),:);
        [fo, G] = fit(lags(inds)', nacf(:,1), ft, fopt);
        data.pro.acfFit = struct('fo',fo,'gof',G, 'fnc',fnc,'res',[lags(inds)' nacf(:,1) - fnc(fo.gamma, fo.tauc, lags(inds)')],'fitCycles',fitCycles);
    end
end

if doPlots
    figure
    clf
    
    for idx = 1:size(acfs,2)
        subplot(1,size(acfs,2),idx)
        inds = find(lags >= 0);
        semilogx(lags(inds), acfs(inds,idx))
        if doFits && idx == 1
            fnc = data.pro.acfFit.fnc;
            fo =  data.pro.acfFit.fo;
            fitCycles = data.pro.acfFit.fitCycles;
            hold on
            if isfield(data.opts,'Vfreq')
                inds = find(lags >= 0 & lags <= fitCycles./data.opts.Vfreq);

                plot(lags(inds), fnc(fo.gamma, fo.tauc, lags(inds)).*acfs(inds(1),idx))
            else
                inds = find(lags >= 0 & lags <= 1);

                plot(lags(inds), fnc(fo.tauc, lags(inds).*acfs(inds(1),idx)))
            end
        end
        xlabel('τ (s)')
        
        if doNorm
            ylabel('NACF')
        else
            ylabel('ACF (m^2)')
        end

        if numel(legCell) >= idx
            legend(legCell{idx})
        end
    end
end

    function fn = get_Fn
        
        if isempty(useField) && isfield(data.opts, 'UseField')
            useField = data.opts.UseField;
            %fprintf('data.opts says use %s for ACF calculation.\n',useField)
        elseif any(strcmp(useField(1), {'x','y'}))
            useField = useField(2:end);
        end
        
        if any(strcmp(direction, {'x', 'y'}))
            if isfield(data.pro, [direction useField])
                fn = useField;
                %fprintf('Going to calculate ACF using field: %s%s\n',direction,useField)
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
                    %fprintf('Going to calculate ACF using field: %s%s\n',d,fn{i})
                elseif isfield(data.pro, [d 'CentresM'])
                    fn{i} = 'CentresM';
                    fprintf('Falling back to calculating ACF using field: %s%s\n',d,fn{i})
                else
                    error('Could not decide a field to use, tried %s and CentresM', useField)
                end
            end
        end
        
    end

end