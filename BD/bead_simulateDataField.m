function data = bead_simulateDataField(opts)
%% data = bead_simulateDataField(opts)
% Wrapper function to call func_simulateDataField and return output as a
% bead data struct

[tracks] = func_simulateDataField(opts);

data = struct('raw', struct(),'opts', struct());
for idx = 1:length(tracks)
    data.raw.xCentresPx(idx,:) = tracks{1,idx}(:,2)';
    data.raw.yCentresPx(idx,:) = tracks{1,idx}(:,3)';
    data.raw.suffixes{idx} = '';
end

data.raw.timeVecMs = tracks{1}(:,1)'*1e3;

data.mPerPx = 1;
data.opts.bdOpts = opts;

data = bead_preProcessCentres(data);