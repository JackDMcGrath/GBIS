function [u v w dwdx dwdy] = mctigue2D(x0,y0,z0,P_G,a,nu,x,y)
% spherical source, forward model based on eq. (52) and (53) by McTigue (1988)
% all parameters are in SI (MKS) units
% u         horizontal (East component) deformation
% v         horizontal (North component) deformation
% w         vertical (Up component) deformation
% dwdx      ground tilt (East component)
% dwdy      ground tilt (North component)
% x0,y0     coordinates of the center of the sphere 
% z0        depth of the center of the sphere (positive downward and
%              defined as distance below the reference surface)
% P_G       dimensionless excess pressure (pressure/shear modulus)
% a         radius of the sphere
% nu        Poisson's ratio
% x,y       benchmark location
% 
% Reference ***************************************************************
% McTigue, D.F. (1987). Elastic Stress and Deformation Near a Finite 
% Spherical Magma Body: Resolution of the Point Source Paradox. J. Geophys.
% Res. 92, 12,931-12,940.
% *************************************************************************
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


% translate the coordinates of the points where the displacement is computed
% in the coordinates system centered in (x0,y0)
xxn = x - x0;
yyn = y - y0;

% radial distance from source center to points where we compute ur and uz
r = sqrt(xxn.^2+yyn.^2);

% dimensionless coordinates
csi = xxn/z0; psi = yyn/z0; rho = r/z0; e = a/z0; 

% constant and expression used in the formulas 
f1 = 1./(rho.^2+1).^1.5;
f2 = 1./(rho.^2+1).^2.5;
c1 = e^3/(7-5*nu);

% displacement (dimensionless) [McTigue (1988), eq. (52) and (53)] 
uzbar = e^3*(1-nu)*f1.*(1 - c1*(0.5*(1+nu)-3.75*(2-nu)./(rho.^2+1)));
urbar = rho.*uzbar;

% displacement (dimensional) 
u = urbar*P_G*z0.*xxn./r;           % East
v = urbar*P_G*z0.*yyn./r;           % North;
w = uzbar*P_G*z0;                   % Up

% GROUND TILT *************************************************************
% see equation (5) in documentation
dwdx = -(1-nu)*P_G*csi*e^3.*f2.*(3-c1*(1.5*(1+nu)-18.75*(2-nu)./(rho.^2+1)));
dwdy = -(1-nu)*P_G*psi*e^3.*f2.*(3-c1*(1.5*(1+nu)-18.75*(2-nu)./(rho.^2+1)));
% *************************************************************************
