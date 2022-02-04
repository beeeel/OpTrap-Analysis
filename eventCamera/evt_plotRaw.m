function evt_plotRaw(cdEvents, varargin)

if nargin == 1
    figure
else
    figure(varargin{1})
end

clf 

subplot(4,2,1)
plot(cdEvents.ts,'.','MarkerSize', 1);
ylabel('Time (μs)')

subplot(4,2,3)
plot(cdEvents.x,'.','MarkerSize', 1);
ylabel('X (px)')

subplot(4,2,5)
plot(cdEvents.y,'.','MarkerSize', 1);
ylabel('Y (px)')

subplot(4,2,7)
if isfield(cdEvents,'gray')
    plot(cdEvents.gray,'.','MarkerSize', 1);
    ylabel('Gray (arb. U)')
elseif isfield(cdEvents,'p')
    plot(cdEvents.p,'.','MarkerSize', 1);
    ylabel('Polarity (+/-)')
end

% Histograms!
subplot(4,2,2)
histogram(cdEvents.ts./60e6)
xlabel('Time (mins)')

ax = subplot(4,2,4);
histogram(cdEvents.x, 'BinLimits', prctile(cdEvents.x,[1 99]));
xlabel('X (px)')
ax.YScale = 'log';

ax = subplot(4,2,6);
histogram(cdEvents.y, 'BinLimits', prctile(cdEvents.y,[1 99]));
xlabel('Y (px)')
ax.YScale = 'log';

subplot(4,2,8)
if isfield(cdEvents,'gray')
    histogram(cdEvents.gray);
    xlabel('Gray (arb. U)')
elseif isfield(cdEvents,'p')
    histogram(cdEvents.p);
    xlabel('Polarity (+/-)')
end