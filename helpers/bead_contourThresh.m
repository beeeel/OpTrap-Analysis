function bead_contourThresh(data)

ths = [];
for idx = 1:length(data.raw.suffixes)
    if startsWith(data.raw.suffixes{idx}, '0')
        str = strsplit(data.raw.suffixes{idx},'th');
        ths(idx) = str2double(str{end});
    end
end

im = max(cat(3,data.Imstack{1}{:,1}),[],3);

figure()
clf
imagesc(im)
axis image
hold on
[c, h] = contour(im, ths,'Color',[1 0 0]);
clabel(c, h)