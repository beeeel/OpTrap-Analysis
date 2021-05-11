function centres = imCentreOfMass2(ims,varargin)
% centres = imCentreOfMass2(ims, [subWidth, method, n])
%% Find 2 centres of mass for each frame in input image array
% Input ims - numeric array of size [X, Y, Z]. Calculates centre of mass
% for each plane z = 1:Z. Second input subWidth (optional) - how many
% pixels from the edge to use for calculating centre of mass. Third input
% method (optional) - 'norm-square' or 'simple'.
%  Returns centres array of size [4 Z] with [XcentreL; YcentreL; XcentreR;
%  YcentreR] for each plane

validateattributes(ims, {'numeric'},{'nonempty','real'},'imCentreOfMass','ims',1)
[imH, imW, nIms] = size(ims);

ims = double(ims);

if nargin <= 2
    method = 'simple';
elseif nargin == 3 || nargin == 4
    method = varargin{2};
else
    disp('Input arguments:')
    disp(varargin)
    error('Too many nargins!')
end
subWidth = imW/2;
if nargin >=2
    subWidth = varargin{1};
end

[X, Y] = meshgrid(cast(1:subWidth,'like',ims),cast(1:imH,'like',ims));

switch method
    case 'simple'
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        Lcentres = squeeze([sum(ims(:,1:subWidth,:).*X,[1 2]); sum(ims(:,1:subWidth,:).*Y,[1 2])]...
            ./sum(ims(:,1:subWidth,:),[1 2]));
        Rcentres = squeeze([sum(ims(:,end-subWidth+1:end,:).*X,[1 2]); sum(ims(:,end-subWidth+1:end,:).*Y,[1 2])]...
            ./sum(ims(:,end-subWidth+1:end,:),[1 2]));
    case 'square'
        % Square every pixel first
        ims = ims.^2;
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        Lcentres = squeeze([sum(ims(:,1:subWidth,:).*X,[1 2]); sum(ims(:,1:subWidth,:).*Y,[1 2])]...
            ./sum(ims(:,1:subWidth,:),[1 2]));
        Rcentres = squeeze([sum(ims(:,end-subWidth+1:end,:).*X,[1 2]); sum(ims(:,end-subWidth+1:end,:).*Y,[1 2])]...
            ./sum(ims(:,end-subWidth+1:end,:),[1 2]));
    case 'norm-square'
        ims = double(ims);
        imsVec = reshape(ims,[],1,nIms);
        imsVecNorm = imsVec - double(mean(ims, [1 2]));
        imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        Lcentres = squeeze([sum(imsMatNorm(:,1:subWidth,:).^2.*X,[1 2]); sum(imsMatNorm(:,1:subWidth,:).^2.*Y,[1 2])]...
            ./sum(imsMatNorm(:,1:subWidth,:).^2,[1 2]));
        Rcentres = squeeze([sum(imsMatNorm(:,end-subWidth+1:end,:).^2.*X,[1 2]); sum(imsMatNorm(:,end-subWidth+1:end,:).^2.*Y,[1 2])]...
            ./sum(imsMatNorm(:,end-subWidth+1:end,:).^2,[1 2]));
    case 'norm-n'
        if nargin == 4
            n = varargin{3};
            validateattributes(n, {'numeric'},{'integer','positive'},'imCentreOfMass','n',2)
        else
            error('When using norm-n, you need another input argument for n')
        end
        ims = double(ims);
        imsVec = reshape(ims,[],1,nIms);
        imsVecNorm = imsVec - double(mean(ims, [1 2]));
        imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        Lcentres = squeeze([sum(imsMatNorm(:,1:subWidth,:).^n.*X,[1 2]); sum(imsMatNorm(:,1:subWidth,:).^n.*Y,[1 2])]...
            ./sum(imsMatNorm(:,1:subWidth,:).^n,[1 2]));
        Rcentres = squeeze([sum(imsMatNorm(:,end-subWidth+1:end,:).^n.*X,[1 2]); sum(imsMatNorm(:,end-subWidth+1:end,:).^n.*Y,[1 2])]...
            ./sum(imsMatNorm(:,end-subWidth+1:end,:).^n,[1 2]));
    case 'dark-norm-square'
        ims = double(ims);
        imsVec = reshape(ims,[],1,nIms);
        imsVecNorm = imsVec - double(mean(ims, [1 2]));
        imsVecNorm(imsVecNorm > 0) = 0;
        imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        Lcentres = squeeze([sum(imsMatNorm(:,1:subWidth,:).^2.*X,[1 2]); sum(imsMatNorm(:,1:subWidth,:).^2.*Y,[1 2])]...
            ./sum(imsMatNorm(:,1:subWidth,:).^2,[1 2]));
        Rcentres = squeeze([sum(imsMatNorm(:,end-subWidth+1:end,:).^2.*X,[1 2]); sum(imsMatNorm(:,end-subWidth+1:end,:).^2.*Y,[1 2])]...
            ./sum(imsMatNorm(:,end-subWidth+1:end,:).^2,[1 2]));
end
Rcentres = Rcentres + [imW - subWidth; 0];
centres = [Lcentres; Rcentres];
