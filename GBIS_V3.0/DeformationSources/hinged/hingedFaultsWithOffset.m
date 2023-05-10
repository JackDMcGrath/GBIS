function [displ] = hingedFaultWithOffset(m,obs,nu, drawgeom)
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
% Last update: 8 August, 2018 - Hinged Dyke
% Last update: 30 July, 2021  - Hinged Fault - Jack McGrath
% Last Update: 3 March, 2022  - XY Location of surface trace - Jack McGrath
% Last Update: 3 October, 2022 - Add fault perp offset (ie different dip of locked) - Jack McGrath

if nargin == 3
    drawgeom=0;
end

% Coordinates of fault trace
X_Trace = m(6);
Y_Trace = m(7);
X0 = 0;
Y0 = 0;

Z1 = m(3);
dip1 = m(4);
Strike = m(5);

X1 = X0+Z1/tand(dip1);
Y1 = Y0;

% Size of X-offset
if length(m) < 14
    staticX=0;
else
    staticX = m(14);
end

% Create vector of Xs and Ys
xs = [X0 X1];
ys = [Y0 Y1];

% Rotate coordinates
xRot = xs*cosd(Strike) - ys*sind(Strike);
yRot = xs*sind(Strike) + ys*cosd(Strike);

X_Fault = xRot(2)+X_Trace;
Y_Fault = -yRot(2)+Y_Trace;

% % Coordinates of upper fault
% X_Fault = m(6); % Center X
% Y_Fault = m(7); % Center Y

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

if drawgeom ==1
figure
drawmodel(Fault1', 'color','r','projection','3D');
hold on
drawmodel(Fault2', 'color','g','projection','3D');
plot(X_Trace,Y_Trace,'g*');
plot3(X_Fault,Y_Fault,-Z1,'r*');
axis equal
xlabel('x axis (m)')
ylabel('y axis (m)')
% view(180-Strike,0)
title(strcat('Original Setup'));
end

% Create vector of Xs and Ys
xs = [X1 X2];
ys = [Y1 Y2];

% Rotate coordinates
xRot = xs*cosd(Strike) - ys*sind(Strike);
yRot = xs*sind(Strike) + ys*cosd(Strike);

% Combine parameters of new rotated Faults
Fault3 = [L1 W1 Z1 dip1 Strike xRot(1)+X_Fault-cosd(Strike)*staticX yRot(1)+Y_Fault+sind(Strike)*staticX ss1 ds1 0];
Fault4 = [L2 W2 Z2 dip2 Strike xRot(2)+X_Fault-cosd(Strike)*staticX -yRot(2)+Y_Fault+sind(Strike)*staticX ss2 ds2 0];

if drawgeom ==1
figure
drawmodel(Fault3', 'color','r','projection','3D');
hold on
drawmodel(Fault4', 'color','g','projection','3D');
plot(X_Trace,Y_Trace,'g*');
plot3(xRot(1)+X_Fault-cosd(Strike)*staticX,yRot(1)+Y_Fault+sind(Strike)*staticX,-Z1,'r*');
axis equal
xlabel('x axis (m)')
ylabel('y axis (m)')
% view(180-Strike,0)
title(strcat('Rotation of ',num2str(Strike),' degrees'));
end

% Initialise displacement matrix
displ = zeros(length(obs(1,:)),3)';

u3 = disloc(Fault3',obs(1:2,:),nu);
u4 = disloc(Fault4',obs(1:2,:),nu);
displ = u3 + u4;
