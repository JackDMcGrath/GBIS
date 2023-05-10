function [displ] = backslip(m,obs,nu)
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
%
% Email: gbis.software@gmail.com
%
% Reference: 
% Bagnardi M. & Hooper A, (2018). 
% Inversion of surface deformation data for rapid estimates of source 
% parameters and uncertainties: A Bayesian approach. Geochemistry, 
% Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
% The function may include third party software.
% =========================================================================
% Last update: 29 October 2021, Jack McGrath

f=[m(1:9);0]; % Parameters for the shear zone (set opening to 0)
Udis = disloc(f,obs,nu); % Slip calculation


Rvel=m(10); % Velocity to right of shear
Lvel=m(11); % Velocity to left of shear

X_center=m(6);
Y_center=m(7);
Strike=m(5);
L=m(1)/2;

% Locations of the ends of the shear 
X1 = X_center-sind(Strike)*L;
Y1 = Y_center-cosd(Strike)*L;

X2 = X_center+sind(Strike)*L;
Y2 = Y_center+cosd(Strike)*L;

[~,s]=point_to_line(obs',[X1,Y1,0],[X2,Y2,0]); %Work out which side of the fault everything is

r=find(s==1); % Observations on right of fault
l=find(s~=1); % Observations on the left of the fault

background=ones(size(obs,2),1);
background(r)=Rvel;
background(l)=Lvel;

U=ones(size(obs,2),3);
U(:,1)=background*sind(Strike);
U(:,2)=background*cosd(Strike);
U(:,3)=0;

% figure
% hold on
% plot(obs(r,1),obs(r,2),'b.')
% plot(obs(l,1),obs(l,2),'k.')
% plot([X1,X2],[Y1,Y2],'r')
% plot(X_center,Y_center,'g*')
% 
% figure
% scatter(obs(:,1),obs(:,2),5,background,'filled');colorbar
% hold on
% plot([X1,X2],[Y1,Y2],'r')
% 
% figure
% scatter(obs(:,1),obs(:,2),5,U(:,1),'filled');colorbar
% hold on
% plot([X1,X2],[Y1,Y2],'r')
% 
% figure
% scatter(obs(:,1),obs(:,2),5,U(:,2),'filled');colorbar
% hold on
% plot([X1,X2],[Y1,Y2],'r')

displ=U'-Udis;