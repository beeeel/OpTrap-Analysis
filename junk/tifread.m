function im = tifread(imfile)
%% Wrapper function to load multiplane TIF files

info = imfinfo(imfile);

im = imread(imfile,'tiff','Info',info);

im(end,end,length(info)) = 0;

for idx = 2:length(info)
    im(:,:,idx) = imread(imfile,'tiff','Info',info,'Index',idx);
end
end