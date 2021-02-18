function centres = imCentreOfMass(ims,varargin)
%% Find centre of mass of each frame in input image array
% Input ims - numeric array of size [X, Y, Z]. Calculates centre of mass
% for each plane z = 1:Z. Second input method (optional) - 'norm-square'
% or 'simple'.
%  Returns centres array of size [2 Z] with [Xcentre; Ycentre] for each
%  plane

validateattributes(ims, {'numeric'},{'nonempty','real'},'imCentreOfMass','ims',1)
[imH, imW, nIms] = size(ims);

ims = double(ims);

if nargin == 1
    method = 'simple';
elseif nargin == 2
    method = varargin{1};
else
    disp('Input arguments:')
    disp(varargin)
    error('Too many nargins!')
end

switch method
    case 'simple'
        [X, Y] = meshgrid(cast(1:imW,'like',ims),cast(1:imH,'like',ims));
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        centres = squeeze([sum(ims.*X,[1 2]) sum(ims.*Y,[1 2])]./sum(ims,[1 2]));
    case 'square'
        % Square every pixel first
        ims = ims.^2;
        [X, Y] = meshgrid(cast(1:imW,'like',ims),cast(1:imH,'like',ims));
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        centres = squeeze([sum(ims.*X,[1 2]) sum(ims.*Y,[1 2])]./sum(ims,[1 2]));
    case 'norm-square'
        ims = double(ims);
        imsVec = reshape(ims,[],1,nIms);
        imsVecNorm = imsVec - double(mean(ims, [1 2]));
        imsMatNorm = reshape(imsVecNorm,imH, imW, nIms);
        [X, Y] = meshgrid(cast(1:imW,'like',ims),cast(1:imH,'like',ims));
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        centres = squeeze([sum(imsMatNorm.^2.*X,[1 2]) sum(imsMatNorm.^2.*Y,[1 2])]./sum(imsMatNorm.^2,[1 2]));
end