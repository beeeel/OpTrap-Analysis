function meta = ParseExperiment_dev(info, meta)
%% Parse optrap experiments by splitting into sections of constant centre coordinate
% Put data into the meta struct for output
%%
Step = 50;
Thresh = 5;

Centres = [info.centres]';
N = hist(Centres);

Diffs = abs(diff([info(1:Step:end).centres]'));

Stds = std([info.centres], 0, 2);

IStart = find(Diffs > Thresh, 1, 'first') * Step;
IEnd = find(Diffs > Thresh, 1, 'last') * Step;

