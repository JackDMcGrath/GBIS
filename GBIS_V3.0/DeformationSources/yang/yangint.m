function [U1,U2,U3] = yangint(x,y,z,z0,theta,a1,b1,a,b,csi,mu,nu,Pdila)
% compute the primitive of the displacement for a prolate ellipsoid
% equation (1)-(8) from Yang et al (JGR, 1988) 
% corrections to some parameters from Newmann et al (JVGR, 2006)
% 
% IN
% x,y,x     coordinates of the point(s) where the displacement is computed [m]
% y0,z0     coordinates of the center of the prolate spheroid (positive downward) [m]
% theta     plunge angle [rad]
% a1,b1     pressure (stress) (output from yangpar.m) [units of P]
% a         semimajor axis [m]
% b         semiminor axis [m]
% c         focus of the prolate spheroid (output from yangpar.m) [m]
% mu        shear modulus [Pa]
% nu        Poisson's ratio
% Pdila     pressure (proportional to double couple forces) [units of P]
% 
%
% OUT
% U1,U2,U3 : displacement in local coordinates [m] - see Figure 3 of Yang et al (1988)
%
% Notes:
% The location of the center of the prolate spheroid is (x0,y0,z0)
%     with x0=0 and y0=0;
% The free surface is z=0;
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


% precalculate parameters that are used often
sint = sin(theta); cost = cos(theta);                                       % y0 = 0;

% new coordinates and parameters from Yang et al (JGR, 1988), p. 4251
% dimensions [m]
csi2 = csi*cost; csi3 = csi*sint;                                           % see Figure 3 of Yang et al (1988)
x1 = x;  x2 = y; x3 = z - z0; xbar3 = z + z0;                               % x2 = y - y0;
y1 = x1; y2 = x2 - csi2; y3 = x3 - csi3; ybar3 = xbar3 + csi3;
r2 = x2*sint - x3*cost; q2 = x2*sint + xbar3*cost;
r3 = x2*cost + x3*sint; q3 = -x2*cost + xbar3*sint;
rbar3 = r3 - csi; qbar3 = q3 + csi;
R1 = sqrt(y1.^2+y2.^2+y3.^2); R2 = sqrt(y1.^2+y2.^2+ybar3.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%C0 = y0*cost + z0*sint;
C0= z0/sint;                                                                % correction base on test by FEM by P. Tizzani IREA-CNR Napoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = (q2*cost + (1+sint)*(R2+qbar3))./(cost*y1 + 1E-15);                  % add 1E-15 to avoid a Divide by Zero warning at the origin

% precalculate parameters that are used often
drbar3 = R1+rbar3; dqbar3 = R2+qbar3; dybar3 = R2+ybar3;
lrbar3 = log(R1+rbar3); lqbar3 = log(R2+qbar3); lybar3 = log(R2+ybar3);
atanb = atan(beta);

% primitive parameters from Yang et al (1988), p. 4252
Astar1 = a1./(R1.*drbar3) + b1.*(lrbar3+(r3+csi)./drbar3);
Astarbar1 = -a1./(R2.*dqbar3) - b1.*(lqbar3+(q3-csi)./dqbar3);

A1 = csi./R1 + lrbar3; Abar1 = csi./R2 - lqbar3;
A2 = R1 - r3.*lrbar3; Abar2 = R2 - q3.*lqbar3;
A3 = csi.*rbar3./R1 + R1; Abar3 = csi.*qbar3./R2 - R2;

Bstar = (a1./R1+2*b1.*A2) + (3-4*nu)*(a1./R2+2*b1.*Abar2);
B = csi*(csi+C0)./R2 - Abar2 - C0*lqbar3;

% the 4 equations below have been changed to improve the fit to internal deformation
Fstar1 = 0; 
Fstar2 = 0; 
F1 = 0; 
F2 = 0; 

f1 = csi*y1./dybar3 + (3/cost^2)*(y1.*sint.*lybar3 - y1.*lqbar3 + 2*q2.*atanb) + 2*y1.*lqbar3 - 4*xbar3.*atanb/cost;
f2 = csi*y2./dybar3 + (3/cost^2)*(q2.*sint.*lqbar3 - q2.*lybar3 + 2*y1*sint.*atanb + cost*(R2-ybar3)) - 2*cost*Abar2 + ...
                    (2/cost)*(xbar3.*lybar3 - q3.*lqbar3);                      % correction after Newmann et al (2006), eq (A-9)
f3 = (1/cost)*(q2.*lqbar3 - q2.*sint.*lybar3 + 2*y1.*atanb) + 2*sint*Abar2 + q3.*lybar3 - csi;


% precalculate coefficients that are used often
cstar = (a*b^2/csi^3)/(16*mu*(1-nu)); cdila = 2*cstar*Pdila;

% displacement components (2) to (7): primitive of equation (1) from Yang et al (1988)
Ustar1 = cstar*(Astar1.*y1 + (3-4*nu)*Astarbar1.*y1 + Fstar1.*y1);          % equation (2) from Yang et al (1988)
% U2star and U3star changed to improve fit to internal deformation
Ustar2 = cstar*(sint*(Astar1.*r2 + (3-4*nu)*Astarbar1.*q2 + Fstar1.*q2) + ...
         cost*(Bstar-Fstar2));                                              % equation (3) from Yang et al (1988)

% The formula used in the script by Fialko and Andy is different from
% equation (4) of Yang et al (1988)
% I use the same to continue to compare the results 2009 07 23
% Ustar3 = cstar*(-cost*(Astarbar1.*r2 + (3-4*nu)*Astarbar1.*q2 - Fstar1.*q2) + ...
%         sint*(Bstar+Fstar2) + 2*cost^2*z.*Astarbar1);           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The equation below is correct - follows equation (4) from Yang et al (1988)
Ustar3 = cstar*(-cost*(Astar1.*r2 + (3-4*nu)*Astarbar1.*q2 - Fstar1.*q2) + ...
         sint*(Bstar+Fstar2));                                              % equation (4) from Yang et al (1988)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Udila1 = cdila*((A1.*y1 + (3-4*nu)*Abar1.*y1 + F1.*y1) - 4*(1-nu)*(1-2*nu)*f1);             % equation (5) from Yang et al (1988)
Udila2 = cdila*(sint*(A1.*r2 + (3-4*nu)*Abar1.*q2 + F1.*q2) - 4*(1-nu)*(1-2*nu)*f2 + ...
         4*(1-nu)*cost*(A2+Abar2) + cost*(A3-(3-4*nu)*Abar3 - F2));                         % equation (6) from Yang et al (1988)
Udila3 = cdila*(cost*(-A1.*r2 + (3-4*nu)*Abar1.*q2 + F1.*q2) + 4*(1-nu)*(1-2*nu)*f3 + ...
         4*(1-nu)*sint*(A2+Abar2) + sint*(A3+(3-4*nu)*Abar3 + F2 - 2*(3-4*nu)*B));          % equation (7) from Yang et al (1988)

    
% displacement: equation (8) from Yang et al (1988) - see Figure 3
U1 = Ustar1 + Udila1;                                                                       % local x component
U2 = Ustar2 + Udila2;                                                                       % local y component
U3 = Ustar3 + Udila3;                                                                       % local z component

