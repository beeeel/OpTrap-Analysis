function datOut = endpoints_norm(accumulated, allplots, tRange, varargin)
%% Pick range of timepoints from results from accumulator_norm
% datOut = endpoints_norm(accumulated, allplots, tRange, [doNorm / deltaT])
% This will pick out all the observations where obsT
% (accumulated{1,n}{2,m}) is within range specified. Specify which
% endpoints you want using allplots. If you supply 4th argument, data is
% normalised by first observation chosen. If 4th argument is numeric, only
% two observations are returned, separated by deltaT = 4th arg.
% 
% datOut is a cell array of size [2, numel(allplots)] where the first row
% is a data array and the second row is the allplot name. The data array
% has size [N, 4], where N is the number of observations. First column is
% observation time, second and third columns are values in each
% dimension and fourth column is observation index (e.g. cell number)

norm = false;
deltaT = NaN;

if nargin == 4
    norm = true;
    if isnumeric(varargin{1})
        deltaT = varargin{1};
    end
end

datOut = cell(2,numel(allplots));

for plt = 1:numel(allplots)
    for dim = 1:2
        
        allDat = [];
        allTs = [];
        alldT = [];
        allcIdx = [];
        cIdx = 1;
        for dayIdx = 1:size(accumulated,2)
            for cellIdx = 1:size(accumulated{1,dayIdx},2)
                clear dat err hack

                t = accumulated{1,dayIdx}{2,cellIdx};
                if isnan(deltaT)
                    % If we're not doing change over time, take all within
                    % range
                    idxs = t > tRange(1) & t < tRange(2);
                else
                    % Otherwise it's more complex - if there's negative
                    % times, take the last one. Otherwise take t(1) and
                    % closest to t(1)+deltaT
                    idxs = false(size(t));
                    if any(t < 0)
                        ind = find(t < 0, 1, 'last');
                    else
                        ind = 1;
                    end
                    idxs(ind) = 1;
                    
                    dt = t - t(ind) - deltaT;
                    [~, ind] = min(abs(dt));
                    idxs(ind) = 1;
                    if sum(idxs) ~= 2
                        error('Some problem with timepoint selection - is 2*deltaT < t(2)-t(1)?')
                    end
                end
                
                if any(idxs)
                    if any(contains(allplots, {'D','α'}))
                        try
                            fps = cat(3,accumulated{1,dayIdx}{1,cellIdx}.fps);
                        catch ME
                            error(ME.message)
                        end
                    end
                    if any(contains(allplots, {'α_τ', 'D_τ', 'α_{τG}', 'D_{τG}', 'D_G'}))
                        try
                            fpstau = cat(3,accumulated{1,dayIdx}{1,cellIdx}.fpstau);
                            fpstauG = cat(3,accumulated{1,dayIdx}{1,cellIdx}.fpstauG);
                            fpsG = cat(3, accumulated{1,dayIdx}{1,cellIdx}.fpsG);
                        catch ME
                            if ~strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch')
                                error(ME.message)
                            end
                        end
                    end
                    if any(contains(allplots, {'α_H', 'D_H', 'α_{τH}', 'D_{τH}', 'α_{τGH}', 'D_{τGH}', 'D_GH'}))
                        try
                            fpsH = cat(3,accumulated{1,dayIdx}{1,cellIdx}.fpsH);
                            fpstauH = cat(3,accumulated{1,dayIdx}{1,cellIdx}.fpstauH);
                            fpstauGH = cat(3,accumulated{1,dayIdx}{1,cellIdx}.fpstauGH);
                            fpsGH = cat(3, accumulated{1,dayIdx}{1,cellIdx}.fpsGH);
                        catch ME
                            if ~strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch')
                                error(ME.message)
                            end
                        end
                    end
                    
                    try % writing robust code and then give up and use try/catch statements instead.
                        switch allplots{plt}
                            case 'α'
                                alphas = reshape(squeeze(fps(1,1,:)), 2, []);
                                dat = alphas(dim,:);
                            case 'α_τ'
                                alphas = reshape(squeeze(fpstau(1,1,:)), 2, []);
                                dat = alphas(dim,:,1);
                            case 'α_{τG}'
                                alphas = reshape(squeeze(fpstauG(1,1,:)), 2, []);
                                dat = alphas(dim,:,1);
                                
                            case 'α_H'
                                alphas = reshape(squeeze(fpsH(1,2,:)), 2, []);
                                dat = alphas(dim,:);
                            case 'α_{τH}'
                                alphas = reshape(squeeze(fpstauH(1,1,:)), 2, []);
                                dat = alphas(dim,:,1);
                            case 'α_{τGH}'
                                alphas = reshape(squeeze(fpstauGH(1,1,:)), 2, []);
                                dat = alphas(dim,:,1);
                                
                            case 'D'
                                Ds = reshape(squeeze(fps(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                            case 'D_τ'
                                Ds = reshape(squeeze(fpstau(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                            case 'D_{τG}'
                                Ds = reshape(squeeze(fpstauG(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                            case 'D_G'
                                Ds = reshape(squeeze(fpsG(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                                
                            case 'D_H'
                                Ds = reshape(squeeze(fpsH(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                            case 'D_{τH}'
                                Ds = reshape(squeeze(fpstauH(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                            case 'D_{τGH}'
                                Ds = reshape(squeeze(fpstauGH(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                            case 'D_{GH}'
                                Ds = reshape(squeeze(fpsGH(2,1,:)), 2, []);
                                dat = Ds(dim,:);
                            
                            case 'M_{1s}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff);
                                dat = Gs(dim,:,1);
                            case 'G_{avg}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff);
                                dat = Gs(dim,:,2);
                            case 'zPos'
                                dat = cat(2,accumulated{1,dayIdx}{1,cellIdx}.zPos);

                            case 'α_{min}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = Gs(dim,:,2);
                            case {'G_{min}', 'βG_0'}
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = Gs(dim,:,1);
                            case '(βG_0)^{-1}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = 1./Gs(dim,:,1);
                            case 'τ_{min}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = Gs(dim,:,3);
                                
                            case 'α_{short}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff3);
                                dat = Gs(dim,:,2);
                            case 'G_{1s}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff3);
                                dat = Gs(dim,:,1);
                            case 'J_{short}'
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff3);
                                dat = Gs(dim,:,3);
                                
                            case 'actual time'
                                dat = t;
                            case 'τ_c'
                                tCs = cat(1, accumulated{1,dayIdx}{1,cellIdx}.tnorm);
                                dat = tCs(:,dim)';
                            case '(τ_c)^{-1}'
                                tCs = cat(1, accumulated{1,dayIdx}{1,cellIdx}.tnorm);
                                dat = 1./tCs(:,dim)';
                            case 'τ_{cH}'
                                tCs = cat(1, accumulated{1,dayIdx}{1,cellIdx}.tnormH);
                                dat = tCs(:,dim)';
                            case 'n_t'
                                dat = cat(2, accumulated{1,dayIdx}{1,cellIdx}.nP);
                            case 'tMax'
                                dat = cat(2, accumulated{1,dayIdx}{1,cellIdx}.tMax);
                                
                            case 'βG_0/τ_c'
                                tCs = cat(1, accumulated{1,dayIdx}{1,cellIdx}.tnorm);
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = log(Gs(dim,:,1))./log(tCs(:,dim))';
                            case 'βG_0/τ_c^2'
                                tCs = cat(1, accumulated{1,dayIdx}{1,cellIdx}.tnorm);
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = log(Gs(dim,:,1))./(log(tCs(:,dim)).^2)';
                            case 'βG_0*τ_c'
                                tCs = cat(1, accumulated{1,dayIdx}{1,cellIdx}.tnorm);
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = (Gs(dim,:,1)).*(tCs(:,dim))';
                            case 'βG_0*τ_c^2'
                                tCs = cat(1, accumulated{1,dayIdx}{1,cellIdx}.tnorm);
                                Gs = cat(2,accumulated{1,dayIdx}{1,cellIdx}.stiff2);
                                dat = log(Gs(dim,:,1)).*(log(tCs(:,dim)).^2)';
                                
                            case 'X'
                                xyz = cat(1,accumulated{1,dayIdx}{1,cellIdx}.XYZ);
                                dat = xyz(:,1);
                            case 'Y'
                                xyz = cat(1,accumulated{1,dayIdx}{1,cellIdx}.XYZ);
                                dat = xyz(:,2);
                            case 'Z'
                                xyz = cat(1,accumulated{1,dayIdx}{1,cellIdx}.XYZ);
                                dat = xyz(:,3);
                                
                            otherwise
                                error('Have you added a new plotName but not how to plot it?')
                        end
                    catch ME
                        if strcmp(ME.identifier, 'MATLAB:badsubscript')
                            warning('Skipped day %s, cell %i, due to subscript error (is there a missing field?)', accumulated{2,dayIdx}, cellIdx)
                            dat = nan(size(idxs));
                        else
                            error(ME.message)
                        end
                    end
                    
                    if length(dat) ~= length(idxs)
                        warning('C''mon you gotta have sensible warning messages. ''a'' is not good enough.')
                    else
                        
                        t = t(idxs);
                        dat = dat(idxs);
                        
                        if norm
                            if isnan(deltaT)
                                dat = dat / dat(1);
                            else
                                dat = dat(2) - dat(1);
                                t = t(2) - t(1);
                            end
                        end
                        
                        
                        allDat = [allDat; dat(:)];
                        allTs = [allTs; t(:)];
                        allcIdx = [allcIdx; repmat(cIdx, numel(dat), 1)];
                        cIdx = cIdx + 1;
                    end
                end
            end
        end

        
        tOut = allTs;
        datOut{1,plt}(:,1) = tOut;
        datOut{1,plt}(:,dim+1) = allDat;
        datOut{2,plt} = allplots{plt};

    end
    datOut{1,plt}(:,end+1) = allcIdx;
end
