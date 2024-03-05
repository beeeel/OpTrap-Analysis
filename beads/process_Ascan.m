function Ascan = process_Ascan(fbase, scanNo, roi, doplots, normalisation, beadPosition, bstop, cropT, tRs, Cref, etaref, zcalibration, Z0V)
dl = dir(sprintf('%sAscan%i*', fbase, scanNo));
fnames = natsort({dl.name}');

if strcmp(normalisation, '0V') && (~exist('Z0V','var') || numel(Z0V) < length(fnames))
    error('Insufficient 0V calibration data')
end

% Check the Ascan threshold
for fIdx = [1 round(length(fnames)/2) length(fnames)]
    data = bead_loadData(fnames{fIdx},true);
    bead_contourThresh(data);
end
    
a = data.opts.beadDiam/2;
deltaZ = logspace(log10(a/2), -4, 1e4);

corPar = faxens_law(a, deltaZ);
corNor = faxens_law(a, deltaZ, false);

Xforce  = nan(length(fnames),1);
Yforce  = nan(length(fnames),1);
Zforce  = nan(length(fnames),1);
Ypos    = nan(length(fnames),1);
Xpos    = nan(length(fnames),1);
Zpos    = nan(length(fnames),1);
Zbar    = nan(length(fnames),1);
ks      = nan(length(fnames),3);
thresh  = nan(length(fnames),1);

First = true;
%%
for fIdx = 1:length(fnames)
    try
        fn = fnames{fIdx};
        data = bead_loadData(fn, true);
        
        [ths, inds] = natsort(data.raw.suffixes(1:end/2));
        try
        % Better way, by regionpropsing the image
        idx = length(ths);
        while idx > 0
        	idx = idx - 1;
            str = strsplit(ths{idx},'th');
            th = str2double(str{2});
            im = max(cat(3,data.Imstack{1}{:,1}),[],3);
            rp = regionprops(im > th);
            if length(rp) > 1 && th ~= 0 && th > mode(im, 'all')
                if idx ~= length(ths)
                    idx = idx + 1;
                    str = strsplit(ths{idx},'th');
                    th = str2double(str{2});
                end
                % Not sorry future me
                data.opts.centresRow = inds(idx);
                thresh(fIdx) = th;
                break
            end 
        end
        catch ME2
            if strcmp(ME2.identifier, 'MATLAB:nonExistentField')
                data.opts.centresRow = inds(idx);
            else
                error('Woopsy!')
            end
        end
            
        data.opts = rmfield(data.opts, 'Vfreq');
        data.opts.bandstop = bstop;
        data.opts.cropT = cropT;
        
        data = bead_preProcessCentres(data);
        
        data = bead_filter_bandstop(data, doplots);
        data = bead_normMSD(data,'doPlots',doplots);
        %     data = bead_PSD(data, 'nblocking', 100,'forceRun',true);
        
        t = data.pro.timeVecMs'*1e-3;
        dc = byteStreamToDouble([data.dirPath sprintf('/I0_th%i.dat',th)]);
        
        bg = byteStreamToDouble([data.dirPath '/I1_th0.dat']);
        dc = dc(1,:);
        
        z = dc(:)./bg(:);
        z = z(data.opts.cropT(1):data.opts.cropT(2));
        switch normalisation
            case 'ref'
                z = z - zbar;
            case 'first'
                if First
                    xbar = mean(data.raw.xCentresPx(data.opts.centresRow,:));
                    ybar = mean(data.raw.yCentresPx(data.opts.centresRow,:));
                    
                    roi = data.opts.roi;
                    
                    First = false;
                end
            case 'none'
                zbar = 0;
                Zbar(fIdx) = mean(z);
            case '0V'
                z = z - Z0V(fIdx);
        end
        
        track = {[t, z]};
        m = msdanalyzer(1, 'Arb.U.', 's','log');
        m = m.addAll(track);
        m = m.computeMSD;
        
        [dydx, tout] = msd_gradientor(m.msd{1}(:,1), m.msd{1}(:,2), 'lsq',5);
        
        
        data = bead_ACF(data,'doFit',true,'forceRun',true, 'doPlots',doplots);
        oCxy = 1./[data.pro.acfFit(1).fo.tauc data.pro.acfFit(2).fo.tauc];
        
        [cT, fP] = msd_cornerator(m, 0, tRs,'doPlot',doplots);
        oCz = 1./cT;
        
        Mxy = 1e-12 * mmsd(data.pro.amsdObj, [0.1 1]);
        Mz = mmsd(m, [1 10]);
        
        kxy = 2*kBT(298) ./ Mxy;
        
        etaxy = kxy ./ (6 * pi * 0.5 * data.opts.beadDiam * oCxy);
        
        etaCorr = etaxy(1) ./ etaref;
        ind = find(corPar < etaCorr,1);
        
        if ~isempty(ind)
            Zpos(fIdx) = deltaZ(ind);
        else
            Zpos(fIdx) = NaN;
        end
        
        if strcmp(zcalibration, 'each') && ~isempty(ind)
            kz = corPar(ind) * oCz * kxy(1) / corNor(ind) / oCxy(1);
            C = sqrt(2 * kBT(298) / Mz / kz);
            %             C = sqrt((Mxy(1)./Mz) * corNor(ind) * oCxy(1) / corPar(ind) / oCz); % numerical error / typo?
        else
            C = Cref;
            kz = 2 * kBT(298) / C^2 / Mz;
        end
        
        
        
        fprintf('Kx = %.3g pN/um, Ky = %.3g pN/um, Kz = %.3g pN/um\nη = %.3g Pa.s\n', 1e6*kxy, 1e6*kz, etaxy(1))
        
        
        if doplots
            figure(4)
            clf
            yyaxis left
            loglog(m.msd{1}(:,1), m.msd{1}(:,2),'k-','LineWidth',2)
            ylabel('MSD (arb.U.)')
            yyaxis right
            semilogx(tout, dydx, 'r--')
            ylabel('Power law α')
            xlabel('\tau (s)')
            
            figure()
            clf
            plot(t, 1e9*z.*C)
            xlabel('Time (s)')
            ylabel('Z (nm)')
        end
        
        if beadPosition
            Ypos(fIdx) = data.opts.XYZ(2) + (mean(data.raw.yCentresPx(data.opts.centresRow,:)) + (1081 - data.opts.roi(2))) * data.mPerPx * 1e6;
            Xpos(fIdx) = data.opts.XYZ(1) + (mean(data.raw.xCentresPx(data.opts.centresRow,:)) + data.opts.roi(1)) * data.mPerPx * 1e6;
            
            Xforce(fIdx) = kxy(1) * (mean(data.raw.xCentresPx(data.opts.centresRow,:)) - xbar - data.opts.roi(1) + roi(1)) * data.mPerPx;
            Yforce(fIdx) = kxy(2) * (mean(data.raw.yCentresPx(data.opts.centresRow,:)) - ybar - data.opts.roi(2) + roi(2)) * data.mPerPx;
        else
            Ypos(fIdx) = data.opts.XYZ(2);
            Xpos(fIdx) = data.opts.XYZ(1);
        end
        ks(fIdx,:) = [kxy, kz];
        
        Zforce(fIdx) = kz * mean(z*C-zbar);
    catch ME
%                     error(ME.message )
    end
    %         close all
end

Ascan.Xforce    = Xforce;
Ascan.Yforce    = Yforce;
Ascan.Zforce    = Zforce;
Ascan.Ypos      = Ypos;
Ascan.Xpos      = Xpos;
Ascan.Zpos      = Zpos;
Ascan.Zbar      = Zbar;
Ascan.ks        = ks;
Ascan.thresh    = thresh;
Ascan.fnames    = fnames;