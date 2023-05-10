function [displ] = splithingedFault(m,obs,nu)
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
% Last update: 28 October, 2021 - Split Hinged Fault - Jack McGrath

Strike = m(1);
X_split = m(2);
Y_split = m(3);

L1 = m(4)/2;
L2 = m(15)/2;
Ld1 = m(5);
Ld2 = m(16);
dip1 = m(8);
dip2 = m(19);

% Locations of centerpoints for surface trace of each split
X1 = X_split-sind(Strike)*L1;
Y1 = Y_split-cosd(Strike)*L1;

X2 = X_split+sind(Strike)*L2;
Y2 = Y_split+cosd(Strike)*L2;

% Locations of centerpoints of each fault at depth
X1 = X1+cosd(Strike)*(Ld1/tand(dip1));
Y1 = Y1-sind(Strike)*(Ld1/tand(dip1));

X2 = X2+cosd(Strike)*(Ld2/tand(dip2));
Y2 = Y2-sind(Strike)*(Ld2/tand(dip2));

% Run these as two hinged faults 
%  [L(m);  W1(m);  Z(m); Dip1(deg); Str(deg); X(m); Y(m);  SS1(m); DS1(m); W2(m); Dip2(deg); SS2(m); DS2(m)]
b1=[m(4);   m(7);  m(5);      m(8);     m(1);   X1;   Y1;    m(9);  m(10); m(11);     m(12);  m(13);  m(14)];
b2=[m(15); m(18); m(16);     m(19);     m(1);   X2;   Y2;   m(20);  m(21); m(22);     m(23);  m(24);  m(25)];

U1 = hingedFaults(b1,obs,nu);
U2 = hingedFaults(b2,obs,nu);

displ = U1 + U2;