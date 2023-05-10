function plothinge(m)
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
%
% Script for only plotting the geometry of hinged okada models
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
% Last update: 30 July, 2021  - Hinged Fault - Jack McGrath

% Coordinates of upper fault
X_Fault = m(6); % Center X
Y_Fault = m(7); % Center Y

% Parameters of upper fault
L1 = m(1); %
W1 = m(2); %
Z1 = m(3); %
dip1 = m(4); %
Strike = m(5);
ss1 = m(8);
ds1 = m(9);
Theta1=0.1;
X1 = 0;
Y1 = 0;

% Parameters of lower fault
L2 = L1;
W2 = m(10); %
Z2 = Z1-W1*sind(dip1);
dip2 = m(11); %
Theta2 = Theta1;
X2 = X1-W1*cosd(dip1);
Y2 = Y1;
ss2 = m(12);
ds2 = m(13);

% Combine parameters into vectors
Fault1 = [L1 W1 Z1 dip1 Theta1 X1+X_Fault Y1+Y_Fault ss1 ds1 0];
Fault2 = [L2 W2 Z2 dip2 Theta2 X2+X_Fault Y2+Y_Fault ss2 ds2 0];

% drawmodel(Fault1', 'color','k','projection','3D'); 
% drawmodel(Fault2', 'color','k','projection','3D');

% Create vector of Xs and Ys
xs = [X1 X2];
ys = [Y1 Y2];

% Rotate coordinates
xRot = xs*cosd(Strike) - ys*sind(Strike);
yRot = xs*sind(Strike) + ys*cosd(Strike);

% Combine parameters of new rotated Faults
Fault3 = [L1 W1 Z1 dip1 Strike xRot(1)+X_Fault yRot(1)+Y_Fault ss1 ds1 0];
Fault4 = [L2 W2 Z2 dip2 Strike xRot(2)+X_Fault -yRot(2)+Y_Fault ss2 ds2 0];

drawmodel(Fault3', 'color','r','projection','no');
drawmodel(Fault4', 'color','g','projection','no');