function U = yangSource(m,obs,nu)

% Function that adapts input and output to use yangdisp.m from dMODELS
% References ***************************************************************
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
a = m(4);
A = m(5);
P_G = m(8);
mu = 1;
theta = m(7); % plunge
phi = m(6);   % strike

x = obs(1,:);
y = obs(2,:);
z = obs(3,:);

% SINGULARITIES ***********************************************************
if theta >= 89.99
    theta = 89.99;                                                          % solution is singular when theta = 90°
end;    
if A >= 0.99
    A = 0.99;
end;  
% *************************************************************************


% DISPLACEMENT ************************************************************
% define parameters used to compute the displacement
b = A*a;                                                                    % semi-minor axis
lambda = 2*mu*nu/(1-2*nu);                                                  % first Lame's elatic modulus
P = P_G*mu;                                                                 % excess pressure
theta = pi*theta/180;                                                       % dip angle in rad
phi = pi*phi/180;                                                           % strike angle in rad

% compute 3D displacements
[u v w] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y,z);
% *************************************************************************

% Combind 3D displacements in U matrix
U = [u;v;w;];