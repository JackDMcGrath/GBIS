function [cmap]=crop_cmap(cmap,clim,center)
% Function to crop divergent cmaps when caxis is not centered

if nargin < 3
    center=0;
end

clim=clim-center;

if abs(clim(2))>abs(clim(1))
    crop=round((256/2)*(1-abs(clim(1)/clim(2))));
    if crop == 0
        crop=1;
    end
    cmap=cmap(crop:end,:);
elseif abs(clim(1))>abs(clim(2))
    crop=round((256/2)*abs(clim(2)/clim(1)));
    cmap=cmap(1:128+crop,:);
end