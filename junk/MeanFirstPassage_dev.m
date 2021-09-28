%% Look at mean first passage time
% First passage time is the time taken to travel more than D distance from
% starting position. Average over all possible starting times.

% Set up a 2D logical operation
% column vector: distance of each observation from chosen starting position
% row vector: chosen distances, Ds.
% Idea 1:
% % result: true if distance < Ds
% % sum trues to find number of timesteps before passing distance D
% Idea 2:
% % result: false if distance < Ds
% % cumsum, logical not, sum to find number of timesteps
% Place into running totals array
% Need sum, sum of squares, number of observations. 
% Calculate mean and std, convert into units of time.

%% First get some data

if ~exist('data','var')
    error('This script needs a data struct to take traces from')
end

x = data.pro.xCentresM(1,:);
y = data.pro.yCentresM(1,:);
t = data.raw.timeVecMs*1e-3;
    
%% Look at this beautiful data
% This took way longer than it should have
figure(2)
clf
xy = 'xy';
rt = {'radial','tangential'};
for idx = 1:2
    subplot(2,2,idx)
%     plot(t, normalize(data.raw.([xy(idx) 'CentresPx'])(1,:)*data.mPerPx*1e6, 'center'))
    plot(t, data.raw.([xy(idx) 'CentresPx'])(1,:)*data.mPerPx*1e6)
    
    xlabel('Time (s)')
    ylabel([xy(idx) ' position (\mum)'])
    title('Raw position')
    
    subplot(2,2,idx+2)
%     plot(t, normalize(eval(xy(idx))*1e6, 'center'))
    plot(t, eval(xy(idx))*1e6)
    
%     plot(t, eval(xy(idx))*1e6-data.raw.([xy(idx) 'CentresPx'])(1,:)*data.mPerPx*1e6)

    xlabel('Time (s)')
    ylabel([rt{idx} ' position (\mum)'])
    title('Angle corrected position')
end

%% Actually do the MFPT calculation

nD = 200;

Ds = logspace(log10(min(diff([x; y])