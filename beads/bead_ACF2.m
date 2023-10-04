function data = bead_ACF2(data, varargin)
%% data = bead_ACF2(data, ...)
% Take data struct and calculate NACF using MSD, returning calculated NACF
% in complete data struct.
%
% Additional parameters in the form of name-value pairs. Possible options:
%   direction       - Direction used for ACF calculation. 'x','y', or 'a'
%   doPlots         - Do the plots if true, else just the calculations
%   centresRow      - Which row of the centres matrix to use. Default all.
%   doNorm          - Normalize (divide) the ACF by the variance of the position
%   useField        - Specify which processed data field to use.
%   forceRun        - Force analysis to run, even if calculation was already done for this data

% Parse the inputs
p = inputParser;

p.addRequired('data',@(x) isa(x,'struct') && isscalar(x) );

p.addParameter('forceRun',false, @(x)islogical(x))
p.addParameter('direction','a',@(x)any(strcmp(x,{'a','all','x','y'})))
p.addParameter('doPlots',true, @(x)islogical(x))
p.addParameter('doFits',isfield(data.opts,'Vfreq'), @(x)islogical(x))
p.addParameter('centresRow',[],@(x)validateattributes(x,{'numeric'},{'positive','integer','<=',numel(data.raw.suffixes)}))
p.addParameter('doNorm',false,@(x)islogical(x))
p.addParameter('useField', [], @(x)any([isfield(data.pro,x),isfield(data.raw,x)]))
p.addParameter('tauMax', 10, @(x) isnumeric(x) && x>0 && isscalar(x))
% p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3))
% p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
% p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(data, varargin{:});

direction = p.Results.direction;
doPlots = p.Results.doPlots;
doFits = p.Results.doFits;
doNorm = p.Results.doNorm;
useField = p.Results.useField;
forceRun = p.Results.forceRun || data.opts.forceRun;
centresRow = p.Results.centresRow;
cropT = data.opts.cropT;
tauMax = p.Results.tauMax;

if isfield(data.pro, 'acf') && ~forceRun
    nacfs = data.pro.nacf(:,2:end);
    tau = data.pro.nacf(:,1);
    if any(strcmp(direction, {'x', 'y'}))
        legCell = {direction};
    else
        legCell = {'X','Y'};
    end
else
    if ~isfield(data.pro,'xCentresM')
        warning('No calibrated centres data found, running preProcessCentres')
        data = bead_preProcessCentres(data);
    end

    if any(strcmp(direction, {'x', 'y'}))
        legCell = {direction};
        X = data.pro.([direction 'CentresM'])';
    else
        legCell = {'X','Y'};
        X = cat(2, data.pro.xCentresM', data.pro.yCentresM');
    end

    N = size(X,1);
    maxDelay = find(data.pro.timeVecMs*1e-3 > tauMax, 1);

    msd = zeros(maxDelay, size(X,2));
    n = zeros(maxDelay, 1);

    for j = 1 : N - 1
        % Square displacement in bulk
        jmax = min(j+maxDelay, N);
        dX = X(j+1:jmax, :) - X(j, :);
        dr2 = dX .* dX; % Each dimension has 1 column

        msd(1:jmax-j,:) = msd(1:jmax-j,:) + dr2;
        n(1:jmax-j) = n(1:jmax-j) + 1;
    end
    msd = msd ./ n;
    msd = msd(n~=0);
    tau = (1:maxDelay) * diff(data.pro.timeVecMs([1 2])*1e-3);
    tau = tau(n~=0);
    nacfs = 1 - 0.5 * msd ./ var(X);

    data.pro.nacf = [tau' nacfs];

    if doFits
        if isfield(data.opts,'Vfreq')
            % Here's the ACF function that you wanna fit to.
            % p = [gamma, tauc]
            fnc = @(gamma, tauc, x) 1 ./ (1+gamma.^2) .* exp(-x / tauc) + gamma.^2./(1+gamma.^2) .* cos(2*pi*data.opts.Vfreq.*x);
            ft = fittype( fnc);
            fopt = fitoptions('method','Nonlin','StartPoint',[0, 0.1],'Lower',[0 0],'Upper',[10, 10]);

            inds = find(tau <= 20./data.opts.Vfreq);
        else
            fnc = @(tauc, x) exp(-x / tauc);
            ft = fittype(fnc);
            fopt = fitoptions('method','Nonlin','StartPoint',[0.1],'Lower',[0],'Upper',[10]);

            inds = 1:length(tau);
        end

        [fo, G] = fit(tau(inds)', nacfs(inds,1), ft, fopt);
        data.pro.acfFit = struct('fo',fo,'gof',G, 'fnc',fnc,'res',[tau(inds)' nacfs(inds,1) - fnc(fo.gamma, fo.tauc, tau(inds)')]);
    end
end

if doPlots
    figure
    clf
    
    for idx = 1:size(nacfs,2)
        subplot(1,size(nacfs,2),idx)
        semilogx(tau, nacfs(:,idx))
        if doFits && idx == 1
            fnc = data.pro.acfFit.fnc;
            fo =  data.pro.acfFit.fo;
            hold on
            if isfield(data.opts,'Vfreq')
                inds = find(tau <= 20./data.opts.Vfreq);

                plot(tau(inds), fnc(fo.gamma, fo.tauc, tau(inds)).*nacfs(1,idx))
            else
                plot(tau, fnc(fo.tauc, tau).*nacfs(1,idx))
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
            fprintf('data.opts says use %s for ACF calculation.\n',useField)
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
