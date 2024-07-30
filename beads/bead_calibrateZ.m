function data = bead_calibrateZ(data, tRs, etaref, doPlots)
%% data = bead_calibrateZ(data, tRs, etaref, doPlots)

a = data.opts.beadDiam/2;
deltaZ = logspace(log10(a), -4, 1e4);

corPar = faxens_law(a, deltaZ);
corNor = faxens_law(a, deltaZ, false);


z = data.raw.zCentresAU;
t = data.pro.timeVecMs*1e-3;

track = {[t', z']};

m = msdanalyzer(1, 'Arb.U.', 's','log');
m = m.addAll(track);
m = m.computeMSD;

data = bead_ACF(data,'doFit',true,'forceRun',true, 'doPlots',doPlots);
oCxy = 1./[data.pro.acfFit(1).fo.tauc data.pro.acfFit(2).fo.tauc];

[cT, ~] = msd_cornerator(m, 0, tRs,'doPlot',doPlots);
oCz = 1./cT;

if ~isfield(data.pro, 'amsdObj')
    data = bead_normMSD(data, 'doPlots',doPlots);
end

Mxy = 1e-12 * mmsd(data.pro.amsdObj, [0.1 1]);
Mz = mmsd(m, [1 10]);

kxy = 2*kBT(data.opts.Temp) ./ Mxy;

etaxy = kxy ./ (6 * pi * 0.5 * data.opts.beadDiam * oCxy);

etaCorr = etaxy(1) ./ etaref;
ind = find(corPar < etaCorr,1);

kz = corPar(ind) * oCz * kxy(1) / corNor(ind) / oCxy(1);
C = sqrt(2 * kBT(298) / Mz / kz);

data.pro.zCentresM = data.raw.zCentresAU * C;
data.pro.viscosity = [etaxy];
data.pro.zCal = C;
data.pro.kz = kz;