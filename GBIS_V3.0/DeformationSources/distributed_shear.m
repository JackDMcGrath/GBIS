function [v]=distributed_shear(dist,shearpar)
%  Script to model fault parallel velocity across a diffuse shear zone, 
%  rather than a single screw dislocation
%  Modified from Tim Wright's Thesis
%  Currently does not support assymetry
%  Inputs:
%   Dist: Distance from center of the shear zone (km)
%   Shearpar:
%     u0: Velocity across the shear zone, with 0 on xmin side (mm/yr)*
%     z:  Depth (km)
%     x0: Center of the shear zone relative to a fault trace (km)
%     w:  Width of the shear zone (km)
% * Distance is measured from the center of the shear zone. At 0km (ie
% center), velocity = 0.5*u0. 0 mm is ALWAYS on the side of -ve distance.
% The only way to change this is to multiply all distances by -1


u0=shearpar(1);
z=shearpar(2);
x0=shearpar(3);
lw=x0-0.5*shearpar(4);
rw=x0+0.5*shearpar(4);

dist=dist-x0;
v=zeros(numel(dist),1);

for i=1:numel(dist)
    v(i)=u0*(pi()*(rw-lw)+2*(dist(i)-lw)*atan2(dist(i)-lw,z)+2*(rw-dist(i))*atan2(dist(i)-rw,z)+z*log(z^2+(dist(i)-rw)^2)-z*log(z^2+(dist(i)-lw)^2))/(2*pi()*(rw-lw));
end

v=v/1000; % Convert displacement to m (for compatibility with the rest of GBIS)
