function centres = volCentreOfMass(ims, pxSz)
% centres = volCentreOfMass(ims)
%% Find centre of mass of each frame in input image array
% Input ims - numeric array of size [X, Y, Z]. Calculates centre of mass
% for whole volume. 
% Returns centres array of size [3 1] with [Xcentre; Ycentre; Zcentre]

validateattributes(ims, {'numeric'},{'nonempty','real'},'imCentreOfMass','ims',1)
[imH, imW, nIms] = size(ims);

if ~exist('pxSz','var')
    pxSz = [1 1 1];
    warning('No pixel calibration provided - assuming isotropic')
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

    [X, Y, Z] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD), cast(1:nIms,'like',imsD));
    % Calculate centre of mass - sum weighted by pixel location, normalized to
    % image total.
    centres = squeeze([sum(imsD.*X.*pxSz(1),"all") sum(imsD.*Y.*pxSz(2),"all") sum(imsD.*Z.*pxSz(3), "all")]./sum(imsD,"all"));
else
    error('This won''t work - the calculation will at least need %gGB of RAM and I''ve not coded a clever way of using less (yet. If you pester me I''ll do it for you',memNeeded)
    blockSize = floor(warnLimit/(8*imW*imH));
    centres = zeros(2,nIms);
    for idxI = 1:blockSize:nIms
        idxF = min(idxI+blockSize, nIms);
        
        imsD = double(ims(:,:,idxI:idxF));
        
        [X, Y] = meshgrid(cast(1:imW,'like',imsD),cast(1:imH,'like',imsD));
        % Calculate centre of mass - sum weighted by pixel location, normalized to
        % image total.
        centres(:,idxI:idxF) = squeeze([sum(imsD.*X,[1 2]) sum(imsD.*Y,[1 2])]./sum(imsD,[1 2]));

    end
end
