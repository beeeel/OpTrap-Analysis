function [MSDnorm, dTs, msd] = normMSD(offset, num_t, timeVec, centres)

% Prepare data to go into msdanalyzer
tracks = cell(length(offset),1);

for idx = 1:length(offset)
    tracks{idx} = [1e-3 .* timeVec(offset(idx) : offset(idx) + num_t - 1)' ...
        1e6 .* centres(offset(idx) : offset(idx) + num_t - 1)'];
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