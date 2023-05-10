function U = pennySource(m,obs,nu)

% Function that adapts input and output to use fialkodisp.m from dMODELS
% References ***************************************************************
% Fialko, Y, Khazan, Y and M. Simons (2001). Deformation due to a 
% pressurized horizontal circular crack in an elastic half-space, with 
% applications to volcano geodesy. Geophys. J. Int., 146, 181–190
% 
% Battaglia, M., Cervelli, P.F., Murray, J.R., 2013a.dMODELS: A MATLAB software 
% package for modeling crustal deformation near active faults and volcanic centers.
% Journal of Volcanology and Geothermal Research 254, 1?4.
% 
% Battaglia, M., Cervelli, P.F., Murray-Muraleda, J.R., 2013b. Modeling 
% crustal deformation?A catalog of deformation models and modeling approaches. 
% U.S. Geological Survey Techniques and Methods, book 13, chap. B1, 96 p.
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
% Last update: 8 August, 2018
%%
% Assign values in m to model parameters
x0 = m(1);
y0 = m(2);
z0 = m(3);
P_G = m(5);
a = m(4);

% Assign obs to x and y
x = obs(1,:);
y = obs(2,:);
z = 0;

% Calculate displacements
[u v w dV] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y,z);

% Combind 3D displacements in U matrix
U = [u;v;w;];