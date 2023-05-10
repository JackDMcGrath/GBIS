function [displ] = splitFault(m,obs,nu,trace)
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
% Last update: 15 October, 2021 - Split Fault - Jack McGrath


% Coordinates of the split at the surface
X_split = m(6);
Y_split = m(7);

% Global fault parameters
L = m(1)/2; % Split assumed to be in center of the fault, so each fault has a length of L/2
Strike = m(5); % Strike of the fault (will be same for both of the sections)
Z = m(3); % Depth of fault tip (ie. locking depth)

% Parameters of one fault
W1 = (m(2)/sind(-m(4)))-(m(3)/sind(-m(4))); % Width of slipping fault, given by distance surface to base of fault, minus distance surface to locking depth
dip1 = m(4); % Dip of this portion of the fault
ss1 = m(8); % Strike-Slip component of this portion of the fault (independent for each section)
ds1 = m(9); % Dip-Slip component of this portion of the fault (independent for each section)
Theta1 = 0.1;
X1= -(m(3)/tand(m(4)));% Midpoint of fault, before rotated off 090 (taking into account surface trace is not the midpoint
Y1=0;
X1=0;
Y1=0;

% Parameters of other fault
W2 = (m(2)/sind(-m(10)))-(m(3)/sind(-m(10))); % Width of slipping fault, given by distance surface to base of fault, minus distance surface to locking depth
dip2 = m(10); % Dip of this portion of the fault
ss2 = m(11);
ds2 = m(12);
Theta2 = Theta1;
X2= -(m(3)/tand(m(10)));% Midpoint of fault, before rotated off 090 (taking into account surface trace is not the midpoint
Y2=0;
X2=0;
Y2=0;

% Combine parameters into vectors
Fault1 = [L W1 Z dip1 Theta1 X1 Y1 ss1 ds1 0];
Fault2 = [L W2 Z dip2 Theta2 X2 Y2 ss2 ds2 0];

% figure
% drawmodel(Fault1', 'color','r','projection','3D');
% hold on
% drawmodel(Fault2', 'color','g','projection','3D');
% axis equal
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% plot(X_split,Y_split,'g*')
% title('Pre-Rotation');

% Create vector of Xs and Ys
xs = [X1 X2];
ys = [Y1 Y2];

% Rotate coordinates
xRot = xs*cosd(Strike) - ys*sind(Strike);
yRot = xs*sind(Strike) + ys*cosd(Strike);

% Combine parameters of new rotated Faults, with co-ordinates so the top of
% the faults meet at the split point
Fault3 = [L W1 Z dip1 Strike xRot(1)+X_split+sind(Strike-180)*L/2 yRot(1)+Y_split+cosd(Strike-180)*L/2 ss1 ds1 0];
Fault4 = [L W2 Z dip2 Strike xRot(2)+X_split-sind(Strike-180)*L/2 yRot(2)+Y_split-cosd(Strike-180)*L/2 ss2 ds2 0];

% Offset fault XY, given the tear is the surface trace not the fault top
Fault3([6 7])=Fault3([6 7])+[(Z/tand(dip1))*cosd(Strike) -(Z/tand(dip1))*sind(Strike)];
Fault4([6 7])=Fault4([6 7])+[(Z/tand(dip2))*cosd(Strike) -(Z/tand(dip2))*sind(Strike)];

% figure
% drawmodel(Fault3', 'color','r','projection','no');
% hold on
% drawmodel(Fault4', 'color','g','projection','no');
% axis equal
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% title(strcat('Rotation of ',num2str(Strike),' degrees'));
% plot(X_split,Y_split,'g*')
% xline(X_split);
% yline(Y_split);
% xlim([-1e5 1e5]),ylim([-5e4 1.5e5])

% Initialise displacement matrix
displ = zeros(length(obs(1,:)),3)';

u3 = disloc(Fault3',obs(1:2,:),nu);
u4 = disloc(Fault4',obs(1:2,:),nu);
displ = u3 + u4;
