function [v]=distributed_shear(dist,shearpar)
%  Script to model fault parallel velocity across a diffuse shear zone, 
%  rather than a single screw dislocation
%  Modified from Tim Wright's Thesis
%  Currently does not support assymetry
%  Inputs:
%   Dist: Distance from center of the shear zone (km)
%   Shearpar:
%     u0: Velocity across the shear zone, with 0 on xmin side (mm/yr)
%     z:  Depth (km)
%     lw: X-Width of shear zone (km, left of 0)
%     rw: X-Width of shear zone(km, right of 0)
%     dx: Fault perpendicular offset of shear zone and fault trace

u0=shearpar(1);
z=shearpar(2);
lw=shearpar(3);
rw=shearpar(4);
dx=shearpar(5);

dist=dist-dx;
v=zeros(numel(dist),1);

for i=1:numel(dist)
    v(i)=u0*(pi()*(rw-lw)+2*(dist(i)-lw)*atan2(dist(i)-lw,z)+2*(rw-dist(i))*atan2(dist(i)-rw,z)+z*log(z^2+(dist(i)-rw)^2)-z*log(z^2+(dist(i)-lw)^2))/(2*pi()*(rw-lw));
end

v=v/1000; % Convert to m