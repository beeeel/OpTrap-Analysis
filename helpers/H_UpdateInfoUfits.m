function info = H_UpdateInfoUfits(info, u_fits)
%% Takes the fits from Unwrap cell and puts them into info.
for fr = 1:length(info)
    info(fr).uMajorAxisLength = u_fits(1,fr);
    info (fr).uMinorAxisLength = u_fits(2,fr);
    info(fr).uOrientation = u_fits(3,fr);
    info(fr).uTaylorParameter = (u_fits(1,fr) - u_fits(2,fr))./(u_fits(1,fr) + u_fits(2,fr));
end
end