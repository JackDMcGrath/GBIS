function [gps]=GPS2FP(gps,bearing,px)
%% GPS2FP
% Reproject E-N GPS velocities into fault parallel - fault perpendicular
% Input: gps: structure containing ve and vn
%        bearing: fault bearing

phi=atand(gps.ve./gps.vn); % Phi = angle of the GPS vector
gpsmag=gps.ve./sind(phi); % Magnitude of GPS vector (from East)
tau=bearing-phi; % tau is angle between GPS and fault
gps.vx=cosd(tau).*gpsmag; % Fault Parallel Velocity
gps.vy=-sind(tau).*gpsmag; % Fault Perp Velocity
gps.vx(isnan(gpsmag))=0; % Set nan values back to 0
gps.vy(isnan(gpsmag))=0; % Set nan values back to 0







end