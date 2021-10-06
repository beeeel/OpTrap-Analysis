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
%
% Both those ideas rely on a large vector operation run for every
% observation. Instead, run a while loop for each initial observation,
% stepping through subsequent observations until dX > Ds. Probably best to
% run through every Ds before moving onto next initial observation.

%% First get some data

if ~exist('data','var')
    error('This script needs a data struct to take traces from')
end

x = data.pro.xCentresM(1,:);
y = data.pro.yCentresM(1,:);
t = data.raw.timeVecMs*1e-3;
dt = mean(diff(t));
    
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

%% MFPT using while loop
% x = sin(2*pi*(0:1e-3:1e1));
tic
nD = 100;
dMin = 1e-9;
dMax = range(x);

Ds = logspace(floor(log10(dMin)), ceil(log10(dMax)), nD);
xFPi = nan(size(x,2)-1,size(Ds,2));
yFPT = zeros(size(x,2)-1,size(Ds,2));
nFPT = 0;
nt = size(x,2);

% Not 100% sure about this: 
% For each initial observation
for idx = 1:size(x,2) - 1
    % Reset continuation flag
    cont = true;
    % Start with first subsequent observation
    jdx = idx + 1;
    % and first distance
    kdx = 1;
    D = Ds(kdx);
    % Keep going through subsequent observations
    while jdx < size(x,2) && cont
        % Calculating displacement
        dx = abs(x(idx) - x(jdx));
        % If the displacement is bigger than the passage distance
        while dx > D && cont
            % Record the first passage indices
            xFPi(idx, kdx) = (jdx - idx);
            % Increment the passage distance
            kdx = kdx + 1;
            cont = kdx <= nD;
            D = Ds(kdx);
            % This will end loop when D > dx or we've done all the passage
            % distances
        end
        % Then it's time to increment the subsequent observation index
        jdx = jdx + 1;
    end
end
        
toc
xMFPT = mean(xFPi.*dt,1,'omitnan');
loglog(Ds, xMFPT)
%% Parallel MFPT using while loop
% x = sin(2*pi*(0:1e-3:1e1));
tic
nD = 100;
dMin = 1e-9;
dMax = range(x);

Ds = logspace(floor(log10(dMin)), ceil(log10(dMax)), nD);
xFPi = nan(nD,size(x,2)-1);
% yFPT = zeros(size(x,2)-1,size(Ds,2));
% nFPT = 0;
% nt = size(x,2);

% Not 100% sure about this: 
% For each initial observation
parfor idx = 1:size(x,2) - 1
    % Reset continuation flag
    cont = true;
    % Initialise first passage index array
    fpi = nan(nD,1);
    % Start with first subsequent observation
    jdx = idx + 1;
    % and first distance
    kdx = 1;
    D = Ds(kdx);
    % Keep going through subsequent observations
    while jdx < size(x,2) && cont
        % Calculating displacement
        dx = abs(x(idx) - x(jdx));
        % If the displacement is bigger than the passage distance
        while dx > D && cont
            % Record the first passage indices
            fpi(kdx) = (jdx - idx);
            % Increment the passage distance
            kdx = kdx + 1;
            cont = kdx <= nD;
            D = Ds(kdx);
            % This will end loop when D > dx or we've done all the passage
            % distances
        end
        % Then it's time to increment the subsequent observation index
        jdx = jdx + 1;
    end
    % Return complete first passage index array
    xFPi(:,idx) = fpi;
end
        
toc
xMFPT = mean(xFPi.*dt,2,'omitnan');
% figure
hold on
loglog(Ds, xMFPT)
xlabel('Distance (m)')
ylabel('MFP time (s)')
% title(['Mean first passage time for ' data.fName])
%% Actually do the MFPT calculation - method 1: parallel/vectorized
tic
nD = 200;
dMin = 1e-9;
dMax = range(x);

Ds = logspace(floor(log10(dMin)), ceil(log10(dMax)), nD);
xFPi = zeros(size(x,2)-1,size(Ds,2));
yFPT = zeros(size(x,2)-1,size(Ds,2));
nFPT = 0;
nt = size(x,2);

% Do the MFPT calculation in parallel
toc
parfor idx = 1:size(x,2)-1
    dx = abs(x(idx+1:nt) - x(idx));
    
    fp = dx' > Ds;
    fps = cumsum(fp);
    fps = ~fps;
    nT = sum(fps);
    
    xFPi(idx, :) = nT * dt;
    
end
xMFPT = mean(xFPi, 2);
toc
