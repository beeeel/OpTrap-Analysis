function bead_imageAndScatter(data)

figure
clf

imagesc(data.ImstackFullFoV{1}{1,1})
axis image
hold on

roi = data.opts.roi;
plot(roi(1) * ones(1,5) + roi(3) * [0 1 1 0 0], roi(2) * ones(1,5) + roi(4) * [0 0 1 1 0], 'r-','LineWidth',2)
plot(data.opts.cCentre(1), data.opts.cCentre(2), 'kx','LineWidth',3)

scatter(data.raw.xCentresPx(1,:)+roi(1), data.raw.yCentresPx(1,:)+roi(2), '.')