function fitout = func_fit_powerspectrum(f, psd, nblock, LfitR, FfitR)
%% Fit PSD, function adapted from Berg-Sorensen tweezercallib

XY = 'XY';
for idx = [1 2]
    P = psd(:,idx); Title = XY(idx);
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        % Bin the powerspectrum
        nbin    =   floor(length(f)/nblock);
        fb = zeros(nbin,1); Pb = zeros(nbin,1); s = zeros(nbin,1);
        for i = 1 : nbin
            fb(i)   = mean(f((i-1)*nblock+1 : i*nblock),'omitnan');
            Pb(i)   = mean(P((i-1)*nblock+1 : i*nblock),'omitnan');
            s(i)    = (1/Pb(i))/sqrt(sum(isfinite(P((i-1)*nblock+1 : i*nblock))));
        end
        fb = fb(:); Pb = Pb(:); s = s(:);
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  Choose data to be fit by a Lorentzian
        ind     = find(fb > LfitR(1) & fb <= LfitR(2));
        xfin    = fb(ind);
        yfin    = Pb(ind);
        sfin    = s(ind);
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  First guess for the fitting parameters FC (corner frequency) and 
        %  D (diffusion coefficient)
        [fc0,D0,sfc,sD,Pfit] = lorentz_analyt(xfin,yfin,nblock);
        
        disp(' ');
        disp(['Lorentzian: fc = ' num2str(fc0,'%5.2f') ' +- ' num2str(sfc,'%5.2f') ', D = ' num2str(D0,'%10.3e') ' +- ' num2str(sD,'%10.3e')]);
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  First guess for F_DIODE (3dB frequency of the photodiode)
        %xxxxxxxxxxxxxxxxx
        %  First, does the diode filter the data, i.e., should we fit f_diode ?

        P_aliasedNyq = sum((D0/(2*pi^2)) ./ ((fNyq + 2*[-10:10]*fNyq).^2 + fc0^2));

        if Pb(length(Pb)) < P_aliasedNyq
            dif         =   Pb(length(Pb)) / P_aliasedNyq;
            f_diode0    =   sqrt(dif * fNyq^2 / (1 - dif));
        else 
            f_diode0   =   2*fNyq;
        end
        
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  Choose data to be fit with the final fit
        
        ind     = find(fb > FfitR(1) & fb <= FfitR(2));
        xfin    = fb(ind);
        yfin    = Pb(ind);
        sfin    = s(ind);
        read_filter_function;
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  Fit    
        check_flag =0;
        weighted_difference = @funn; % Creating a function handle better than defining a string (see matlab help on feval)
        if (want_alpha == 1)
            alpha0 = 0.3;
            a0 = sqrt(1/alpha0^2 - 1);      % Substitute a0 for alpha to ensure that alpha lies between 0 and 1
        else 
            alpha0 = []; 
            a0 =[];    
        end;
        parameters0 = [fc0 D0 f_diode0 a0];             %Initial fitting parameters
        scal_fit = ones(1,length(parameters0));         %Scaled fitting parameters
        disp(' ');
        if (length(parameters0) == 2)
            disp(['Parameters0: fc = ' num2str(parameters0(1),'%5.1f') ', D = ' num2str(parameters0(2),'%10.3e')]);
        elseif (length(parameters0) == 3)
            disp(['Parameters0: fc = ' num2str(parameters0(1),'%5.1f') ', D = ' num2str(parameters0(2),'%10.3e') ',  fdiode = ' num2str(parameters0(3),'%10.3e')]);
        elseif (length(parameters0) == 4)
            disp(['Parameters0: fc = ' num2str(parameters0(1),'%5.1f') ', D = ' num2str(parameters0(2),'%10.3e') ',  fdiode = ' num2str(parameters0(3),'%10.3e') ',    alpha = ' num2str(1/sqrt(1+parameters0(4)^2),'%5.3f')]);
        end
        [scal_fit,RESNORM,RESIDUAL,JACOBIAN] = fit_nonl(weighted_difference,scal_fit,tolx,50,parameters0,xfin,yfin,sfin,check_flag);
        scal_fit = abs(scal_fit);               % The function P_theor is symmetric in alpha and fdiode 
        parameters = scal_fit.*parameters0;
        if (length(parameters) == 2)
            disp(['Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.3e')]);
        elseif (length(parameters) == 3)
            disp(['Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.3e') ',  fdiode = ' num2str(parameters(3),'%10.3e')]);
        else 
            disp(['Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.3e') ',  fdiode = ' num2str(parameters(3),'%10.3e') ',    alpha = ' num2str(parameters(4),'%5.3f')]);
        end
        disp(['chi^2 = ',num2str(RESNORM,'%6.2f')]);
        nfree = length(yfin) - length(parameters0);
        bbac = 1. - gammainc(RESNORM/2.,nfree/2.);     %Calculate backing of fit
        chi2 = RESNORM/nfree;
        disp(['chi^2 per degree of freedom = ',num2str(chi2,'%6.2f') ', n_{free} = ' num2str(nfree,'%6.2f')]);
        disp(['backing = ',num2str(bbac,'%10.3e')]);
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  If nblock is small, add extra minimization criteria
        if nblock < 200,
            parameters0 = parameters;
            scal_fit = ones(1,length(parameters0));
            sumfun =  inline('sum(((1./P_theor(scal_fit,parameters0,xfin,check_flag) - 1./yfin).^2) ./ (sfin.^2) + 2*log(P_theor(scal_fit,parameters0,xfin,check_flag)))',...
                'scal_fit','parameters0','xfin','yfin','sfin','check_flag');       
            [scal_fit,RESNORM,EXITFLAG,OUTPUT] = fminsearch(sumfun,scal_fit,...
                optimset('Display','iter','MaxIter',100,'MaxFunEvals',1000,'LargeScale','off'),parameters0,xfin,yfin,sfin,check_flag);
            scal_fit = abs(scal_fit); % The function P_theor is symmetric in alpha and fdiode 
            parameters = scal_fit.*parameters0;
            disp(' ');
            disp(['Second fit. Parameters: fc = ' num2str(parameters(1),'%5.1f') ', D = ' num2str(parameters(2),'%10.1e')]);
            if length(parameters0) > 2, disp([',  fdiode = ' num2str(parameters(3),'%10.1e')]); end;
            if length(parameters0) > 3, disp(['            alpha = ' num2str(1/sqrt(1+parameters(4)^2),'%5.2f')]); end;
            disp(['Minimized sum = ',num2str(RESNORM,'%6.1f')]);
            %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %   Calculate the real chi^2 and backing:
            [scal_fit,RESNORM,RESIDUAL,JACOBIAN] = fit_nonl(weighted_difference,scal_fit,tolx,0,parameters0,xfin,yfin,sfin,check_flag);
            nfree = length(yfin) - length(parameters0);
            bbac = 1. - gammainc(RESNORM/2.,nfree/2.);
            chi2 = RESNORM/nfree;
            disp(['chi^2 per degree of freedom = ',num2str(chi2,'%6.2f') ', n_{free} = ' num2str(nfree,'%6.2f')]);
            disp(['backing = ',num2str(bbac,'%8.0f')]);
        end;
        
        if (length(parameters0) > 3)
            parameters(4) = 1/sqrt(1+parameters(4)^2);
            for i= 1: length(parameters0)-1
                JACOBIAN(:,i)=JACOBIAN(:,i)/parameters0(i);     %Rescaling Jacobian back to unscaled parameters    
            end
            JACOBIAN(:,4) = JACOBIAN(:,4)/(-1/sqrt(1/alpha0-1)*1/sqrt(1/parameters(4)-1)*1/parameters(4)^2);     %Rescaling Jacobian back to unscaled parameters    
        else
            for i= 1: length(parameters0)
                JACOBIAN(:,i)=JACOBIAN(:,i)/parameters0(i);     %Rescaling Jacobian back to unscaled parameters    
            end
        end 
        eval(['scal_fit' Title ' = scal_fit;']);
        eval(['parameters' Title ' = parameters;']);
        eval(['parameters0' Title ' = parameters0;']);
        eval(['bac' Title ' = bbac;']);
        eval(['chi2' Title ' = chi2;']);
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  Errors on parameters
        
        CURVATURE   = JACOBIAN' * JACOBIAN;
        cov = COVARIANCES(CURVATURE);       %This way of inverting the CURVATURE matrix prevents a nearly singlular matrix in many cases
        SIGMA_PAR   = [];
        sigma_par = sigmapar(CURVATURE,parameters);  
        eval(['sigma_par' Title ' = sigma_par;']);   
        eval(['cov' Title ' = cov;']);   
        disp('                                                                                              ')
        disp(['cov(fc,D)            = ',num2str(cov(1)/sqrt(parameters(1)*parameters(2)),'%8.3f')]);
        if  length(parameters0) > 2,
            disp(['cov(fc,fdiode)       = ' num2str(cov(2)/sqrt(parameters(1)*parameters(3)),'%8.3f')]);
            disp(['cov(D,fdiode)        = ' num2str(cov(3)/sqrt(parameters(2)*parameters(3)),'%8.3f')]);
        end;
        if length(parameters0) > 3, 
            disp(['cov(fc,alpha)        = ' num2str(cov(4)/sqrt(parameters(1)*parameters(4)),'%8.3f')]); 
            disp(['cov(D,alpha)         = ' num2str(cov(5)/sqrt(parameters(2)*parameters(4)),'%8.3f')]);
            disp(['cov(fdiode,alpha)    = ' num2str(cov(6)/sqrt(parameters(3)*parameters(4)),'%8.3f')]);
        end;
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        %  Plot the powerspectrum    
        
        if idx == 1, 
            color_now = color_x;
            figure(11); clf; subplot('position',[0.1 .6 .35 .3]); set(gcf,'Numbertitle','off','Name','Consistency of fit'); hold on; h = title('Final Fit, X','FontWeight','bold'); ...
                set(h,'FontUnits','Normalized','FontSize',0.07); plot_P_cos(f,P,scal_fit,parameters0,color_now);
            figure(11); subplot('position',[0.1 .15 .35 .3]); hold on; h = title('Final Fit, X','FontWeight','bold');ind = find(f > Ffit_start & f < Ffit_end);...
            set(h,'FontUnits','Normalized','FontSize',0.07);
            plot_data_div_fit(f(ind),P(ind),scal_fit,parameters0,color_now);
            figure(10); clf; set(gcf,'Numbertitle','off','Name','Final Fit, X'); hold on; h = title('Final Fit, X','FontWeight','bold');ind = find(f > Plot_start & f < Plot_end);...
                plot_fit(f(ind),P(ind),scal_fit,parameters0,color_now,1);
        else 
            color_now = color_y;
            figure(11); subplot('position',[0.6 .6 .35 .3]); hold on; h = title('Final Fit, Y','FontWeight','bold');plot_P_cos(f,P,scal_fit,parameters0,color_now);
            set(h,'FontUnits','normalized','FontSize',0.07);
            figure(11); subplot('position',[0.6 .15 .35 .3]); hold on; h = title('Final Fit, Y','FontWeight','bold');ind = find(f > Ffit_start & f < Ffit_end);...
            set(h,'FontUnits','normalized','FontSize',0.07);
            plot_data_div_fit(f(ind),P(ind),scal_fit,parameters0,color_now);
            figure(13); clf; set(gcf,'Numbertitle','off','Name','Final Fit, Y'); hold on; h = title('Final Fit, Y','FontWeight','bold');ind = find(f > Plot_start & f < Plot_end);...
            plot_fit(f(ind),P(ind),scal_fit,parameters0,color_now,2);
        end;
        
        %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
    end; %(for ixy = [1 2])