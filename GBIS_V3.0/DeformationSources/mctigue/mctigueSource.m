function U = mctigueSource(m,obs,nu)

% Function that adapts input and output to use mctigue2D.m from dMODELS
% References ***************************************************************
% McTigue, D.F. (1987). Elastic Stress and Deformation Near a Finite 
% Spherical Magma Body: Resolution of the Point Source Paradox. J. Geophys.
% Res. 92, 12,931-12,940.
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
% Markov chain Monte Carlo algorithm incorporating the Metropolis alghoritm
% (e.g., Mosegaard & Tarantola, JGR,(1995).
%
% by Marco Bagnardi and Andrew Hooper (COMET, University of Leeds)
% Email: M.Bagnardi@leeds.ac.uk
% Reference: TBA (Bagnardi and Hooper, in prep.)
%
% The function uses third party software.
% =========================================================================
% Last update: 10/05/2017
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

% Calculate displacements
[u v w dwdx dwdy] = mctigue2D(x0,y0,z0,P_G,a,nu,x,y);

% Combind 3D displacements in U matrix
U = [u;v;w;];