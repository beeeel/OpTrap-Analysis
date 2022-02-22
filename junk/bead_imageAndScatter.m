function bead_imageAndScatter(data, fh)

if ~exist('fh','var')
    fh = figure;
else
    figure(fh.Number)
end
clf
fn = strsplit(data.dirPath,'/');
fh.Name = sprintf('Image and scatter for %s',fn{end});

% To get a black and white image, first choose display range by percentile
% and then give imagesc an [m,n,3] matrix.
im = data.ImstackFullFoV{1}{1,1};
lims = prctile(im, [0.1 99.9],'all');
im = im - lims(1);
im(im > lims(2)-lims(1)) = lims(2)-lims(1);
im = double(im)./double(lims(2));

imagesc(im.*ones(1,1,3))
axis image
hold on

roi = data.opts.roi;
plot(roi(1) * ones(1,5) + roi(3) * [0 1 1 0 0], roi(2) * ones(1,5) + roi(4) * [0 0 1 1 0], 'r-','LineWidth',2)

if isfield(data.opts, 'cCentre')
    plot(data.opts.cCentre(1), data.opts.cCentre(2), 'kx','LineWidth',3)
end

scatter(data.raw.xCentresPx(1,:)+roi(1), data.raw.yCentresPx(1,:)+roi(2), [], data.raw.timeVecMs*1e-3, '.')