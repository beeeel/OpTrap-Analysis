function doubleVec = byteStreamToDouble(filePath)
%% Load a double vector from a file containing Java bytestream
%doubleVec = byteStreamToDouble(filePath)
%
% Tested on files saved from fast_acq.bsh on B09 trapping rig
% Due to endianness, the stream needs to be reversed twice to recover
% original data

FID = fopen(filePath);
stream = uint8(fread(FID)');
fclose(FID);
streamReverse = stream(end:-1:1);
doubleVecReverse = typecast(streamReverse,'double');
doubleVec = doubleVecReverse(end:-1:1);

% dname = 'E:/Will/beads/2020_11_09_2um/';
% contents = dir(dname);
% out = struct;
% for file = contents(~[contents.isdir])'
%     FID = fopen([dname file.name]);
%     a = uint8(fread(FID)');
%     fclose(FID)
%     a = a(end:-1:1);
%     c = typecast(a,'double');
%     c = c(end:-1:1);
%     out.(file.name(1:end-4)) = c;
% end