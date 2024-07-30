function data = bead_PSDZ(data, varargin)
%% data = bead_PSDZ(data, ...)
% Take data struct and calculate PSD for processed or raw Z data, returning
% calculated PSD in complete data struct.
%
% Additional parameters in the form of name-value pairs. Possible options:
%   doPlots         - Do the plots if true, else just the calculations
%   doFit           - Do fitting as per Berg-Sørensen 2004, plot fits.
%   forceRun        - Force analysis to run, even if calculation was already done for this data
%   nblocking       - Block averaging after fourier transform - improves SNR, default 20
%   plotAx          - Plot on a pair of given axes
%   zeroPadding     - Zero padding for FFT

% Excess options unused:
%   direction       - Direction used for ACF calculation. 'x','y', or 'a'
%   useField        - Specify which processed data field to use.
%   doNorm          - Normalize (divide) the ACF by the variance of the position
%   centresRow      - Which row of the centres matrix to use. Default all.
%   useRaw          - Use raw data instead of processed (you probably don't want this)

% Parse the inputs
p = inputParser;

p.addRequired('data',@(x) isa(x,'struct') && isscalar(x) );

p.addParameter('forceRun',false, @(x)islogical(x))
p.addParameter('doPlots',true, @(x)islogical(x))
p.addParameter('doFit',true, @(x)islogical(x))
p.addParameter('legCell',{'X','Y'}, @(x)iscell(x))
p.addParameter('nBlocking',20, @(x) isscalar(x) && round(x)==x && x>0)
p.addParameter('plotAx',[], @(x)isa(x,'matlab.graphics.axis.Axes'))
p.addParameter('zeroPadding',[], @(x) isscalar(x) && round(x)==x && x>0)

% p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3))
% p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
% p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(data, varargin{:});

doPlots = p.Results.doPlots;
forceRun = p.Results.forceRun || data.opts.forceRun;
legCell = p.Results.legCell;
nB = p.Results.nBlocking;
doFit = p.Results.doFit;
zp = p.Results.zeroPadding;

centresRow = min(data.opts.centresRow, size(data.pro.zCentresM,1));

% Don't run unless necessary
if isfield(data.pro, 'psdz') && ~forceRun
    psds = data.pro.psdz(:,2:end);
    freq = data.pro.psdz(:,1);
else
    tracks = data.pro.zCentresM(centresRow,:);
    t = data.pro.timeVecMs./1e3;

    % PSD method based on Berg-Sørensen 2004.
    %     [psds, freq] = dft1(tracks', t);
    %     psds = abs(psds(end/2:end,:)).^2 ./t(end); % Take one-sided FT to calculate PSD = FT^2 / tmax
    %     freq = freq(end/2:end);
    
    delta_t =   t(2)-t(1);
    fNyq    =   1 / (2 * delta_t);
    
    T       =   max(t);
    freq       =   ((1 : length(tracks)) / T)';
    
    if zp == 1
        ZP = length(t);
    else
        ZP = 2^nextpow2(max(zp,length(t)));
    end
    FT      =   delta_t*fft(tracks',ZP);
    P       =   FT .* conj(FT) / T;
    
    ind     =   find(freq <= fNyq); % only to the Nyquist f
    ind     = ind(2:end);
    freq       =   freq(ind);
    psds       =   P(ind,:);

    if nB > 1
        % Blocking from Berg-Sørensen 2004 section IV.
        inds = 1:nB:size(psds,1);
        if inds(end) ~= size(psds,1)
            inds(end) = size(psds,1);
        end
        psdsb = zeros(length(inds)-1, size(psds,2));
        wb = zeros(length(inds)-1,1);

        for idx = 1:length(inds)-1
            psdsb(idx,:) = mean(psds(inds(idx):inds(idx+1)-1,:),1);
            wb(idx) = mean(freq(inds(idx):inds(idx+1)-1,:),1);
        end
        psds = psdsb;
        freq = wb;

        clear psdsb wb
    end

    data.pro.psdz = [freq psds];

    if doFit
        tmax = diff(data.pro.timeVecMs([1 end])) / 1e3;
        spq = @(p, q) sum(freq.^(2*p) .* psds.^q, 1);

        fc1 = sqrt( ( spq(0,1) .* spq(2,2) - spq(1,1) .* spq(1,2) ) ...
            ./ ( spq(1,1) .* spq(0,2) - spq(0,1) .* spq(1,2) ) );
        D1 = ( (2*pi^2) / tmax ) * ( spq(0,2) .* spq(2,2) - spq(1,2).^2 ) ...
            ./ ( spq(1,1) .* spq(0,2) - spq(0,1) .* spq(1,2) );

        data.pro.psdFitsz = [fc1; D1];
    end
end

if doPlots
    makeNewAxis = isempty(p.Results.plotAx) || size(psds,2) ~= numel(p.Results.plotAx);
    if makeNewAxis
        figure
        clf
        for idx = 1:size(psds,2)
            ax(idx) = subplot(1,size(psds,2),idx);
        end
    else
        ax = p.Results.plotAx;
    end
    
    set(ax,'XScale','log','YScale','log')
    
    for idx = 1:size(psds,2)
        hold(ax(idx),'on')
        loglog(ax(idx),freq, abs(psds(:,idx)))
        xlabel(ax(idx),'Frequency (Hz)')
        
        ylabel(ax(idx),'PSD (m^2Hz)')
        
        if doFit
            % From Berg-Sørensen 2004 eq 10.
            psd = @(f, D, fc) D ./ (2 * pi^2 * (fc.^2 + f.^2));

            tmax = diff(data.pro.timeVecMs([1 end])) / 1e3;
            spq = @(p, q) sum(freq.^(2*p) .* psds.^q, 1);

            fc1 = sqrt( ( spq(0,1) .* spq(2,2) - spq(1,1) .* spq(1,2) ) ...
                ./ ( spq(1,1) .* spq(0,2) - spq(0,1) .* spq(1,2) ) );
            D1 = ( (2*pi^2) / tmax ) * ( spq(0,2) .* spq(2,2) - spq(1,2).^2 ) ...
                ./ ( spq(1,1) .* spq(0,2) - spq(0,1) .* spq(1,2) );
            plot(ax(idx), freq, psd(freq, D1(idx), fc1(idx)))
        end

        if numel(legCell) >= idx
            legend(legCell{idx})
        end
    end
end

if isfield(data.opts,'Vfreq') && isfield(data.opts,'Vpp')
        dI = 1; % direction index
        
        frange = [-0.5 0.5] + data.opts.Vfreq;
        inds = find(freq > frange(1) & freq < frange(2));
        [P,ind] = max(psds(inds,dI));
        idx = inds(ind);
        if data.opts.Vpp ~= 0 && doPlots
            plot(ax(dI), freq(idx), P, 'kx')
        end
        data.pro.psdAmpz = [freq(idx) psds(idx,:)];
    end

end

function [X, w] = dft1(x,t)
% From Berg-Sørensen 2004.
if ~exist('t','var')
    t = (0:length(x)-1)';
end
dt = diff(t([1 2]));
N = length(t);
j = (1:N)';
X = zeros(size(x));
k = (-N/2+1:N/2)';
w = k ./ (N.*dt);

% for kdx = 1:length(k)
%     X(kdx,:) = dt .* sum(exp(1i*2*pi*j*k(kdx)/N) .*x, 1);
% end
for jdx = 1:length(j)
    X = X + dt .* exp(1i*2*pi*j(jdx)*k/N) .* x(jdx,:);
end
end


% 
% function [X, w] = dft2(x,t)
% if ~exist('t','var')
%     t = (0:length(x)-1)';
% end
% dt = diff(t([1 2]));
% N = length(t);
% j = (1:N)';
% X = zeros(size(x));
% k = -N/2+1:N/2;
% w = k' ./ (N.*dt);
% 
% for kdx = 1:length(k)
%     X(kdx,:) = dt .* sum(exp(1i*2*pi*j*k(kdx)/N) .*x, 1);
% end
% % for jdx = 1:length(j)
% %     X = X + dt .* sum(exp(1i*2*pi*j(jdx)*k/N) .* x(jdx), 1);
% % end
% end
