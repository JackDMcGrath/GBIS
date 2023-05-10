function [Ux Uy Uz] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y,z)
% compute the 3D displacement due to a pressurized ellipsoid  
%
% IN
% a         semimajor axis [m]
% b         semiminor axis [m]
% lambda    Lame's constant [Pa]
% mu        shear modulus [Pa]
% nu        Poisson's ratio 
% P         excess pressure (stress intensity on the surface) [pressure units]
% x,y,x     coordinates of the point(s) where the displacement is computed [m]
% x0,y0,z0  coordinates of the center of the prolate spheroid (positive downward) [m]
% theta     plunge angle [rad]
% phi       trend angle [rad]
%
% OUT
% Ux,Uy,Uz  displacement 
%
% Note ********************************************************************
% compute the displacement due to a pressurized ellipsoid 
% using the finite prolate spheroid model by from Yang et al (JGR,1988)
% and corrections to the model by Newmann et al (JVGR, 2006).
% The equations by Yang et al (1988) and Newmann et al (2006) are valid for a
% vertical prolate spheroid only. There is and additional typo at pg 4251 in 
% Yang et al (1988), not reported in Newmann et al. (2006), that gives an error 
% when the spheroid is tilted (plunge different from 90°):
%           C0 = y0*cos(theta) + z0*sin(theta)
% The correct equation is 
%           C0 = z0/sin(theta)
% This error has been corrected in this script.
% *************************************************************************
%==========================================================================
% USGS Software Disclaimer 
% The software and related documentation were developed by the U.S. 
% Geological Survey (USGS) for use by the USGS in fulfilling its mission. 
% The software can be used, copied, modified, and distributed without any 
% fee or cost. Use of appropriate credit is requested. 
%
% The USGS provides no warranty, expressed or implied, as to the correctness 
% of the furnished software or the suitability for any purpose. The software 
% has been tested, but as with any complex software, there could be undetected 
% errors. Users who find errors are requested to report them to the USGS. 
% The USGS has limited resources to assist non-USGS users; however, we make 
% an attempt to fix reported problems and help whenever possible. 
%==========================================================================


% testing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all; clc;
% a = 1000; b = 0.99*a;
% lambda = 1; mu = lambda; nu = 0.25; P = 0.01;
% theta = pi*89.99/180; phi = 0;
% x = linspace(0,2E4,7);
% y = linspace(0,1E4,7);
% x0 = 0; y0 = 0; z0 = 5E3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the parameters for the spheroid model
[a1,b1,c,Pdila,Pstar] = yangpar(a,b,lambda,mu,nu,P);

% translate the coordinates of the points where the displacement is computed
% in the coordinates systen centered in (x0,0)
xxn = x - x0;
yyn = y - y0;

% rotate the coordinate system to be coherent with the model coordinate
% system of Figure 3 (Yang et al., 1988)
xxp  = cos(phi)*xxn - sin(phi)*yyn;
yyp  = sin(phi)*xxn + cos(phi)*yyn;


% compute displacement for a prolate ellipsoid at csi = c
[U1p,U2p,U3p] = yangint(xxp,yyp,z,z0,theta,a1,b1,a,b,c,mu,nu,Pdila);
% compute displacement for a prolate ellipsoid at csi = -c
[U1m,U2m,U3m] = yangint(xxp,yyp,z,z0,theta,a1,b1,a,b,-c,mu,nu,Pdila);
Upx = -U1p-U1m;
Upy = -U2p-U2m;
Upz  =  U3p+U3m;
% rotate horizontal displacement back (strike)
Ux = cos(phi)*Upx + sin(phi)*Upy;
Uy = -sin(phi)*Upx + cos(phi)*Upy;
Uz = Upz;



