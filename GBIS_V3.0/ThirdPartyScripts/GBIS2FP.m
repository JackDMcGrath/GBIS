function [modGPS]=GBIS2FP(modGPS,bearing)

mod.ve=modGPS(1,:)';
mod.vn=modGPS(2,:)';
[mod]=GPS2FP(mod,bearing,[]);
modGPS(1,:)=mod.vx;
modGPS(2,:)=mod.vy;

end