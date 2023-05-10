function [displ] = hingedDikes(m,obs,nu)
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
%
% Create a 'hinged' dike by connecting two okada models together via a
% shared top-bottom edge
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
% Last update: 8 August, 2018 - Hinged Dyke

% Coordinates of upper dike
X_dike = m(1); %
Y_dike = m(2); %

% Parameters of upper dike
L1 = m(3); %
W1 = m(4); %
Z1 = m(5); %
Phi1 = m(6); %
Theta1 = 0.01;
X1 = 0;
Y1 = 0;
ss1 = 0;
ds1 = 0;
op1 = m(7); %

% Parameters of lower sill
L2 = L1;
W2 = m(8); %
Z2 = Z1-W1*sind(Phi1);
Phi2 = m(9); %
Theta2 = Theta1;
X2 = X1-W1*cosd(Phi1);
Y2 = Y1;
ss2 = 0;
ds2 = 0;
op2 = m(10); %

% Dikes strike
Strike = m(11); %

% Combine parameters into vectors
Dike1 = [L1 W1 Z1 Phi1 Theta1 X1+X_dike Y1+Y_dike ss1 ds1 op1];
Dike2 = [L2 W2 Z2 Phi2 Theta2 X2+X_dike Y2+Y_dike ss2 ds2 op2];

% drawmodel(Dike1', 'color','k','projection','3D'); 
% drawmodel(Dike2', 'color','k','projection','3D');

% Create vector of Xs and Ys
xs = [X1 X2];
ys = [Y1 Y2];

% Rotate coordinates
xRot = xs*cosd(Strike) - ys*sind(Strike);
yRot = xs*sind(Strike) + ys*cosd(Strike);

% Combine parameters of new rotated dikes
Dike3 = [L1 W1 Z1 Phi1 Strike xRot(1)+X_dike yRot(1)+Y_dike ss1 ds1 op1];
Dike4 = [L2 W2 Z2 Phi2 Strike xRot(2)+X_dike -yRot(2)+Y_dike ss2 ds2 op2];

% drawmodel(Dike3', 'color','r','projection','3D'); 
% drawmodel(Dike4', 'color','r','projection','3D');
% axis equal
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% title(strcat('Rotation of',num2str(Strike),' degrees'));

% Initialise displacement matrix
displ = zeros(length(obs(1,:)),3)';

u3 = disloc(Dike3',obs(1:2,:),nu);
u4 = disloc(Dike4',obs(1:2,:),nu);
displ = u3 + u4;
