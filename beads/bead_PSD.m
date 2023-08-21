function data = bead_PSD(data, varargin)
%% data = bead_PSD(data, ...)
% Take data struct and calculate PSD for processed or raw data, returning
% calculated PSD in complete data struct. If ACF has not been calculated,
% do that first.
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
p.addParameter('doPlots',true, @(x)islogical(x))
p.addParameter('legCell',{'X','Y'}, @(x)iscell(x))

% p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3))
% p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
% p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(data, varargin{:});

doPlots = p.Results.doPlots;
forceRun = p.Results.forceRun || data.opts.forceRun;
legCell = p.Results.legCell;

if isfield(data.pro, 'psd') && ~forceRun
    psds = data.pro.psd(:,2:end);
    freq = data.pro.psd(:,1);
else

    if ~isfield(data.pro, 'acf')
        data = bead_ACF(data,'doPlots',false);
    end

    acfs = data.pro.acf(:,2:end);
    lags = data.pro.acf(:,1);

    psds = zeros((size(lags,1)+1)/2,size(acfs,2));
    for idx = 1:size(acfs,2)
        [freq, psds(:,idx)] = fft_scaled(lags, acfs(:,idx), false);
    end

    data.pro.psd = [freq' psds];
end

if doPlots
    figure
    clf
    
    for idx = 1:size(psds,2)
        subplot(1,size(psds,2),idx)
        loglog(freq, abs(psds(:,idx)))
        xlabel('Frequency (Hz)')
        
        ylabel('PSD (m^2/Hz)')
        
        if numel(legCell) >= idx
            legend(legCell{idx})
        end
    end
end

end