function data = bead_normMSD_polyfit(data, field, offset, varargin)

% [MSDnorm, dTs, msd]
% offset, num_t, timeVec, centres)


timeVec = data.raw.timeVecMs;
centres = data.raw.(field);

if nargin == 3
    num_t = length(centres);
elseif nargin == 4
    num_t = varargin{1};
else
    error('Wrong number of input arguments')
end
% offset = [2e4 num_t+1+2e4 2.19e5 num_t+1+2.19e5 1.7505e6 1.7505e6+num_t+1];%1:1e5:1.8e6+1;

% [MSDnorm, dTs, msd] = normMSD(offset, num_t, timeVec, xCentresHP);


% Dimension order for polynomial fitting
dims = [1, 3, 2];

%%
% Prepare data to go into msdanalyzer
tracks = cell(length(offset),1);

for idx = 1:length(offset)
    % Crop data, then fit polynomial of order set by processing script
    centres = centres(offset(idx) : offset(idx) + num_t - 1);
    [~, centres, ~] = func_thermal_rm(1:length(centres), ...
        permute(centres, dims), data.opts.pOrder, 1, length(centres));
    centres = ipermute(centres, dims);
    
    tracks{idx} = [1e-3 .* timeVec(offset(idx) : offset(idx) + num_t - 1)' ...
        1e6 .* centres'];
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

data.pro.([field(1) 'MSDnorm']) = [dTs MSDnorm];
data.pro.([field(1) 'msdObj']) = msd;
%%
fh = figure;
fh.Name = data.fName;
clf
hold on

loglog(dTs, MSDnorm,'LineWidth',2)

fh.Children.XAxis.Scale = 'log';
fh.Children.YAxis.Scale = 'log';

% xlim([1e-3 10])

xlabel('Delay (s)')
ylabel(['Normalized MSD in ' field(1)])
title({'Normalized mean square displacement', ['From t = ' num2str(diff(tracks{1}([1, end], 1))) 's of observations']})
% legend('20 - 40s','40 - 60s', '3m00s - 3m20s','3m20s - 3m40s','32m00s - 32m20s', '32m20s - 32m40s','Location','best')
