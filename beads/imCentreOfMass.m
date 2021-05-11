function centres = imCentreOfMass(ims,varargin)
% centres = imCentreOfMass(ims, [method, n])
%% Find centre of mass of each frame in input image array
% Input ims - numeric array of size [X, Y, Z]. Calculates centre of mass
% for each plane z = 1:Z. Second input method (optional) - 'norm-square'
% or 'simple'.
%  Returns centres array of size [2 Z] with [Xcentre; Ycentre] for each
%  plane

validateattributes(ims, {'numeric'},{'nonempty','real'},'imCentreOfMass','ims',1)
[imH, imW, nIms] = size(ims);

if nargin == 1
        method = 'simple';
    elseif nargin == 2 || nargin == 3
        method = varargin{1};
    else
        disp('Input arguments:')
        disp(varargin)
        error('Too many nargins!')
end

%% Gotta check it's gonna be ok to do it vectorized
memNeeded = 8 * numel(ims); % Number of bytes to store ims as double

warnLimit = 2e9;    % Warning at 2GB seems sensible?
maxLimit = 4e9;     % 4GB max size seems sensible?
if memNeeded < maxLimit
    % Do it vectorized
    if memNeeded > warnLimit
        warning('Might need lots of memory, sorry if I crash!')
    end
    
    imsD = double(ims);

    switch method
        case 'simple'
            [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
            % Calculate centre of mass - sum weighted by pixel location, normalized to
            % image total.
            centres = squeeze([sum(imsD.*X,[1 2]) sum(imsD.*Y,[1 2])]./sum(imsD,[1 2]));
        case 'square'
            % Square every pixel first
            imsD = imsD.^2;
            [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
            % Calculate centre of mass - sum weighted by pixel location, normalized to
            % image total.
            centres = squeeze([sum(imsD.*X,[1 2]) sum(imsD.*Y,[1 2])]./sum(imsD,[1 2]));
        case 'norm-square'
            imsD = double(imsD);
            imsVec = reshape(imsD,[],1,nIms);
            imsVecNorm = imsVec - double(mean(imsD, [1 2]));
            imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
            [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
            % Calculate centre of mass - sum weighted by pixel location, normalized to
            % image total.
            centres = squeeze([sum(imsMatNorm.^2.*X,[1 2]) sum(imsMatNorm.^2.*Y,[1 2])]./sum(imsMatNorm.^2,[1 2]));
        case 'norm-n'
            if nargin == 3
                n = varargin{2};
                validateattributes(n, {'numeric'},{'integer','positive'},'imCentreOfMass','n',2)
            else
                error('When using norm-n, you need another input argument for n')
            end
            imsD = double(imsD);
            imsVec = reshape(imsD,[],1,nIms);
            imsVecNorm = imsVec - double(mean(imsD, [1 2]));
            imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
            [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
            % Calculate centre of mass - sum weighted by pixel location, normalized to
            % image total.
            centres = squeeze([sum(imsMatNorm.^n.*X,[1 2]) sum(imsMatNorm.^n.*Y,[1 2])]./sum(imsMatNorm.^n,[1 2]));
        case 'dark-norm-square'
            imsD = double(imsD);
            imsVec = reshape(imsD,[],1,nIms);
            imsVecNorm = imsVec - double(mean(imsD, [1 2]));
            imsVecNorm(imsVecNorm > 0) = 0;
            imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
            [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
            % Calculate centre of mass - sum weighted by pixel location, normalized to
            % image total.
            centres = squeeze([sum(imsMatNorm.^2.*X,[1 2]) sum(imsMatNorm.^2.*Y,[1 2])]./sum(imsMatNorm.^2,[1 2]));
    end
else
    blockSize = floor(warnLimit/(8*imW*imH));
    allcentres = zeros(2,nIms);
    for idxI = 1:blockSize:nIms
        idxF = min(idxI+blockSize, nIms);
        
        imsD = double(ims(:,:,idxI:idxF));
        
        switch method
            case 'simple'
                [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
                % Calculate centre of mass - sum weighted by pixel location, normalized to
                % image total.
                centres = squeeze([sum(imsD.*X,[1 2]) sum(imsD.*Y,[1 2])]./sum(imsD,[1 2]));
            case 'square'
                % Square every pixel first
                imsD = imsD.^2;
                [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
                % Calculate centre of mass - sum weighted by pixel location, normalized to
                % image total.
                centres = squeeze([sum(imsD.*X,[1 2]) sum(imsD.*Y,[1 2])]./sum(imsD,[1 2]));
            case 'norm-square'
                imsD = double(imsD);
                imsVec = reshape(imsD,[],1,nIms);
                imsVecNorm = imsVec - double(mean(imsD, [1 2]));
                imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
                [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
                % Calculate centre of mass - sum weighted by pixel location, normalized to
                % image total.
                centres = squeeze([sum(imsMatNorm.^2.*X,[1 2]) sum(imsMatNorm.^2.*Y,[1 2])]./sum(imsMatNorm.^2,[1 2]));
            case 'norm-n'
                if nargin == 3
                    n = varargin{2};
                    validateattributes(n, {'numeric'},{'integer','positive'},'imCentreOfMass','n',2)
                else
                    error('When using norm-n, you need another input argument for n')
                end
                imsD = double(imsD);
                imsVec = reshape(imsD,[],1,nIms);
                imsVecNorm = imsVec - double(mean(imsD, [1 2]));
                imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
                [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
                % Calculate centre of mass - sum weighted by pixel location, normalized to
                % image total.
                centres = squeeze([sum(imsMatNorm.^n.*X,[1 2]) sum(imsMatNorm.^n.*Y,[1 2])]./sum(imsMatNorm.^n,[1 2]));
            case 'dark-norm-square'
                imsD = double(imsD);
                imsVec = reshape(imsD,[],1,nIms);
                imsVecNorm = imsVec - double(mean(imsD, [1 2]));
                imsVecNorm(imsVecNorm > 0) = 0;
                imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
                [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
                % Calculate centre of mass - sum weighted by pixel location, normalized to
                % image total.
                centres = squeeze([sum(imsMatNorm.^2.*X,[1 2]) sum(imsMatNorm.^2.*Y,[1 2])]./sum(imsMatNorm.^2,[1 2]));
        end
        allcentres(:,idxI:idxF) = centres;
    end
    centres = allcentres;
end
