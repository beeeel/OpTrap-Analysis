function data = bead_Ddist_and_ngp(data)
%% data = bead_Ddist_and_ngp(data)
% Calculate increment distribution and non-gaussian parameter for bead data

% NGP - see Bursac 2005 Nat.Mat. or Weeks 2002 PRL.
%  <z^4> / (3 <z^2> ^2 ) - 1
ngp = @(z, rho) sum(rho.*(z.^4),1) ./ (3 .* sum(rho.*(z.^2),1) .^2) - 1;

if isfield(data.pro,'amsdObj')
    m = data.pro.amsdObj;
else
    if isfield(data.opts, 'UseField')
        fld = data.opts.UseField;
    else
        fld = 'CentresM';
    end
    if isfield(data.opts, 'CentresRow')
        cRow = data.opt.CentresRow;
    else
        cRow = 1;
    end
    m = msdanalyzer(1, 'um','s','log');
    t = {[data.raw.timeVecMs'*1e-3, data.pro.(['x' fld])(cRow,:)']; 
        [data.raw.timeVecMs'*1e-3, data.pro.(['y' fld])(cRow,:)']};
    m = m.addAll(t);
end
    
m = m.computeDdist;

dt = round(m.Ddist{1,1}(:,1),1, 'significant');
alpha = zeros(size(dt,1), 2);

% Calculate ngp
for dim = 1:2
    rho = m.Ddist{dim,2}(1:floor(end/2),:);
    rho = rho./sum(rho);
    edg = m.Ddist{dim,2}(ceil(end/2):end,1);
    z   = edg(1:end-1) + diff(edg(1:2))/2;
    alpha(:,dim) = ngp(z, rho)';
end

% Alpha plot
figure(300)
h = semilogx(dt, alpha(:,1), 'x-');
hold on
semilogx(dt, alpha(:,2), 'o-', 'Color',h.Color)
xlabel('τ (s)')
ylabel('α')

% ρ(z) plot
fh = figure;
clf
y = normpdf(z);
y = y./sum(y);
% semilogy(z,y, 'k-','LineWidth',2)
hold on
for idx = 1:size(rho,2)
    semilogy(z,rho(:,idx) ./ y,'o-')
end
set(gca,'YScale','log')
fh.Name = data.fName;
xlim([-5 5])
ylim([1e-1 100])

data.pro.ngp = [dt alpha];
