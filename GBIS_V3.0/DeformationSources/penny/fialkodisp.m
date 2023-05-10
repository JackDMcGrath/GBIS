function [u v w dV] = fialkodisp(x0,y0,z0,P_G,a,nu,x,y,z)
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
% *************************************************************************
% Fialko, Y, Khazan, Y and M. Simons (2001). Deformation due to a 
% pressurized horizontal circular crack in an elastic half-space, with 
% applications to volcano geodesy. Geophys. J. Int., 146, 181–190
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

% General parameters ******************************************************
eps = 1E-8;                                                                 % relative accuracy 
rd = a;                                                                     % avoid conflict with line 49 
h = z0/rd;                                                                  % dimensionless source depth                
% *************************************************************************

% Coordinates transformation **********************************************
x = (x-x0)/rd; y = (y-y0)/rd; z = (z-z0)/rd;                                % translate and scale locations
r = sqrt(x.^2+y.^2);                                                        % compute radial distance 
% *************************************************************************

% solve for PHI and PSI, Fialko et al. (2001), eq. (26) *******************
[csi1 w1] = gauleg(eps,10,41);                                                 
[csi2 w2] = gauleg(10,60,41);                                                  
csi = cat(2,csi1,csi2);  wcsi = cat(2,w1,w2);                               % ascissas and weights for Gauss-Legendre quadrature 
if size(csi,1)==1, csi = csi'; end;                                         % check that csi is a column vectors
    
[phi psi t wt] = psi_phi(h);
PHI = sin(csi*t)*(wt'.*phi);                                                % Gauss-Legendre quadrature
PSI = (sin(csi*t)./(csi*t) - cos(csi*t))*(wt'.*psi);                        % Gauss-Legendre quadrature
% *************************************************************************

% compute A and B, Fialko et al. (2001), eq. (24) *************************
% NOTE there is an error in eq (24), (1-exp(-2*a)) must be replaced by exp(-a) 
a = csi*h;
A = exp(-a).*(a.*PSI+(1+a).*PHI);
B = exp(-a).*((1-a).*PSI-a.*PHI);
% *************************************************************************

% compute Uz and Ur, Fialko et al. (2001), eq. (12) and (13) **************
% NOTE there are two errors in eq (12) and (13)
% (1) 2*Uz and 2*Ur must be replaced by Uz and Ur
% (2) dcsi/sinh(csi*h) must be replaced by dcsi
Uz = zeros(size(r));                                                        % pre-allocate variable
Ur = zeros(size(r));                                                        % pre-allocate variable
for i=1:length(r)
    J0 = besselj(0,r(i)*csi);
    Uzi = J0.*(((1-2*nu)*B - csi*(z+h).*A).*sinh(csi*(z+h)) + ... 
                        (2*(1-nu)*A - csi*(z+h).*B).*cosh(csi*(z+h)));
    Uz(i) = wcsi*Uzi;               
    J1 = besselj(1,r(i)*csi);
    Uri = J1.*(((1-2*nu)*A + csi*(z+h).*B).*sinh(csi*(z+h)) + ... 
                        (2*(1-nu)*B + csi*(z+h).*A).*cosh(csi*(z+h)));
    Ur(i) = wcsi*Uri;    
end                   
% *************************************************************************

% Deformation components **************************************************
u = rd*P_G*Ur.*x./r;
v = rd*P_G*Ur.*y./r;
w = -rd*P_G*Uz;
% *************************************************************************

% Volume change ***********************************************************
dV = -4*pi*(1-nu)*P_G*rd^3*(t*(wt'.*phi));
% *************************************************************************