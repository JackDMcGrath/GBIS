function [u v w dV dwdx dwdy eea gamma1 gamma2] = fialko01(x0,y0,z0,P_G,a,nu,x,y,z)
% 3D Green's function for sill-like source (Fialko et al., 2001)
% all parameters are in SI (MKS) units
%
% INPUT
% x0,y0     coordinates of the center of the sphere 
% z0        depth of the center of the sill (positive downward and
%              defined as distance below the reference surface)
% P_G       dimensionless excess pressure (pressure/shear modulus)
% a         radius of the sphere
% nu        Poisson's ratio
% x,y       benchmark location
% z         depth within the crust (z=0 is the free surface)
% 
% OUTPUT
% u         horizontal (East component) deformation
% v         horizontal (North component) deformation
% w         vertical (Up component) deformation
% dV        volume change
% dwdx      ground tilt (East component)
% dwdy      ground tilt (North component)
% eea       areal strain
% gamma1    shear strain
% gamma2    shear strain
%
% Reference ***************************************************************
% Fialko, Y, Khazan, Y and M. Simons (2001). Deformation due to a 
% pressurized horizontal circular crack in an elastic half-space, with 
% applications to volcano geodesy. Geophys. J. Int., 146, 181–190
% *************************************************************************
%
% Note ********************************************************************
% compute the displacement due to a pressurized sill-like source 
% using the finite penny-crack model by Fialko et al (GJI,2001)
% There are two typos in the published paper
% eq (12) and (13)
%  (1) 2*Uz and 2*Ur must be replaced by Uz and Ur
%  (2) dcsi/sinh(csi*h) must be replaced by dcsi
% eq (24)
%      (1-exp(-2*a)) must be replaced by exp(-a) 
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



% DISPLACEMENT ************************************************************
% compute 3D displacements
[u v w dV] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y,z);
% *************************************************************************


% TILT ********************************************************************
h = 0.001*abs(max(x)-min(x));                                              % finite difference step

% East comonent
[~, ~, wp] = fialkodisp(x0,y0,z0,P_G,a,nu,x+h,y,z);
[~, ~, wm] = fialkodisp(x0,y0,z0,P_G,a,nu,x-h,y,z);
dwdx = 0.5*(wp - wm)/h;

% North component
[~, ~, wp] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y+h,z);
[~, ~, wm] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y-h,z);
dwdy = 0.5*(wp - wm)/h;
% *************************************************************************


% STRAIN ******************************************************************
% Displacement gradient tensor
[up , ~, ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x+h,y,z);
[um , ~, ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x-h,y,z);
dudx = 0.5*(up - um)/h;

[up , ~, ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y+h,z);
[um , ~, ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y-h,z);
dudy = 0.5*(up - um)/h;

[~, vp , ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x+h,y,z);
[~, vm , ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x-h,y,z);
dvdx = 0.5*(vp - vm)/h;

[~, vp , ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y+h,z);
[~, vm , ~] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y-h,z);
dvdy = 0.5*(vp - vm)/h;

% Strains
eea = dudx + dvdy;                                                          % areal strain
gamma1 = dudx - dvdy;                                                       % shear strain
gamma2 = dudy + dvdx;                                                       % shear strain
% *************************************************************************
