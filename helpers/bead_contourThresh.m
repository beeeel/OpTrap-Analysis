function bead_contourThresh(data)

ths = [];
for idx = 1:length(data.raw.suffixes)
    if startsWith(data.raw.suffixes{idx}, '0')
        str = strsplit(data.raw.suffixes{idx},'th');
        ths(idx) = str2double(str{end});
    end
end
figure()
clf
imagesc(data.Imstack{1}{1,1})
hold on
[c, h] = contour(data.Imstack{1}{1,1}, ths,'Color',[1 0 0]);
clabel(c, h)