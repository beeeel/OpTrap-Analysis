function [DErrors, Dvals] = SimulateUnwrapFitting(Ia, Dvals, varargin)
%% Calculate the minimum deformation threshold to be detectable above noise
% This needs input parsing to be added - prority is "bin edges" from
% histcounts


% First find at standard deviations for the first 30 frames
tol = 0.15;
n_repeats = 50;
n_frames = 100;
n_theta = 360;
n_reps = 5;

Par = ParseInputs(Ia, Dvals, varargin{:});
n_sims = length(Par.Dvals);

idxa = Ia > (1-tol)*median(Ia,2) & Ia < (1+tol)*median(Ia,2); % Indexes for values included in fitting

fiteqn = @(a, b, phi, x) sqrt((a .* cos(0.0175.*x + phi)).^2 + (b .* sin(0.0175.*x + phi)).^2 );


theta = linspace(single(1),single(360),n_theta);
sim_fits = zeros(3,n_sims,n_repeats);
R = sqrt(mean(Ia(idxa(:,:,1:n_frames)),'all'));
orientn = - pi/2 + pi*rand(1,n_sims,n_repeats);
rng('shuffle')
devs = zeros(n_frames,1);

for frame = 1:n_frames
    devs(frame) = std(Ia(:,idxa(:,:,frame),frame));
end
% Assume this deviation is the noise on a circle
noise = mean(devs);
%%
% Create some fictitious data using the radial equation and noise
% comparable to the noise found above, for a range of D values
tic
fprintf('\n')
for repeat = 1:n_repeats
    for sim = 1:n_sims
        D = Par.Dvals(sim);
        sim_data = fiteqn(R,R*(1-D)./(1+D),orientn(1,sim,repeat),...
            linspace(1,360,n_theta)) + normrnd(0,noise,1,n_theta);
        th_fit = repmat(theta,1,n_reps) + ...
            reshape((0:n_reps-1).*ones(n_theta,1),1,[]);
        fitobj = fit(double(th_fit'),repmat(sim_data',n_reps,1),fiteqn,...
            'Lower', [0 0 -pi/2], 'Upper', [inf inf pi/2], 'Start', [10, 10, 3]);
        sim_fits(:,sim,repeat) = [fitobj.a, fitobj.b, fitobj.phi];
    end
    MyProgressBar(repeat / n_repeats);
end
% Fix the equivalent fitting equations problem - if b>a
sim_fits(3,:) = sim_fits(3,:) + (sim_fits(2,:) > sim_fits(1,:)) * pi/2;
sim_fits(3,:) = mod(sim_fits(3,:)+pi/2, pi) - pi/2;
sim_fits(1:2,:) = [max(sim_fits(1:2,:)); min(sim_fits(1:2,:))];
sim_Dval = abs((sim_fits(1,:,:) - sim_fits(2,:,:)) ./ (sim_fits(1,:,:) + sim_fits(2,:,:)));
DErrors = mean(abs(sim_Dval-Par.Dvals) ./ Par.Dvals,3);

fprintf('%s\n',repmat(' ',1,104))
fprintf('Simulated data in %gs\n',toc)


function Results = ParseInputs(Ia, Dvals, varargin)

Fname = 'SimulateUnwrapFitting';
p = inputParser();

addRequired(p,'Ia',@(x)validateattributes(x,{'numeric'},...
    {'positive','nonempty','nrows',1,'3d'},Fname,'Ia'))
addRequired(p,'Dvals',@(x)validateattributes(x,{'numeric'},...
    {'nonempty','nonnegative','row'},Fname,'Dvals'))
addParameter(p,'NoiseFrames',30,@(x)validateattributes(x,{'numeric'},...
    {'positive','nonempty','<',NFrames}))



parse(p, Ia, Dvals, varargin{:})
Results = p.Results;
end
end
