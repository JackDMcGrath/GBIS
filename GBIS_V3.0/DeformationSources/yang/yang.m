function [u v w dwdx dwdy eea gamma1 gamma2] = yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z)
% 3D Green's function for a spheroidal source 
% all parameters are in SI (MKS) units
%
% OUTPUT
% u         horizontal (East component) deformation
% v         horizontal (North component) deformation
% w         vertical (Up component) deformation
% dwdx      ground tilt (East component)
% dwdy      ground tilt (North component)
% eea       areal strain
% gamma1    shear strain
% gamma2    shear strain
%
% SOURCE PARAMETERS
% a         semimajor axis
% A         geometric aspect ratio [dimensionless]
% P_G       dimennsionless excess pressure (pressure/shear modulus) 
% x0,y0     surface coordinates of the center of the prolate spheroid
% z0        depth of the center of the sphere (positive downward and
%              defined as distance below the reference surface)
% theta     plunge (dip) angle [deg] [90 = vertical spheroid]
% phi       trend (strike) angle [deg] [0 = aligned to North]
%
% CRUST PARAMETERS
% mu        shear modulus
% nu        Poisson's ratio 
%
% BENCHMARKS
% x,y       benchmark location
% z         depth within the crust (z=0 is the free surface)
%
% Reference ***************************************************************
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


% % TILT ********************************************************************
% h = 0.001*abs(max(x)-min(x));                                              % finite difference step
% 
% % East comonent
% [~, ~, wp] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x+h,y,z);
% [~, ~, wm] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x-h,y,z);
% dwdx = 0.5*(wp - wm)/h;
% 
% % North component
% [~, ~, wp] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y+h,z);
% [~, ~, wm] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y-h,z);
% dwdy = 0.5*(wp - wm)/h;
% % *************************************************************************
% 
% 
% % STRAIN ******************************************************************
% % Displacement gradient tensor
% [up , ~, ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x+h,y,z);
% [um , ~, ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x-h,y,z);
% dudx = 0.5*(up - um)/h;
% 
% [up , ~, ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y+h,z);
% [um , ~, ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y-h,z);
% dudy = 0.5*(up - um)/h;
% 
% [~, vp , ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x+h,y,z);
% [~, vm , ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x-h,y,z);
% dvdx = 0.5*(vp - vm)/h;
% 
% [~, vp , ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y+h,z);
% [~, vm , ~] = yangdisp(x0,y0,z0,a,b,lambda,mu,nu,P,theta,phi,x,y-h,z);
% dvdy = 0.5*(vp - vm)/h;
% 
% % Strains
% eea = dudx + dvdy;                                                          % areal strain
% gamma1 = dudx - dvdy;                                                       % shear strain
% gamma2 = dudy + dvdx;                                                       % shear strain
% % *************************************************************************
