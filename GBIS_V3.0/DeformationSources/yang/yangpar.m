function [a1,b1,c,Pdila,Pstar] = yangpar(a,b,lambda,mu,nu,P)
% compute the parameters for the spheroid model
% formulas from [1] Yang et al (JGR,1988)
% corrections from [2] Newmann et al (JVGR, 2006), Appendix
%
% IN
% a         semimajor axis [m]
% b         semiminor axis [m]
% lambda    Lame's constant [Pa]
% mu        shear modulus [Pa]
% nu        Poisson's ratio 
% P         excess pressure (stress intensity on the surface) [pressure units]
%
% OUT
% a1, b1    pressure (stress) [units of P] from [1]
% c         prolate ellipsoid focus [m]
% Pdila     pressure (proportional to double couple forces) [units of P] from [1]
% Pstar     pressure [units of P]
%
% Notes:
% [-]   : dimensionless
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


epsn = 1E-10;

c = sqrt(a^2-b^2);                                  % prolate ellipsoid focus [m]

a2 = a^2; a3 = a^3; b2 = b^2;
c2 = c^2; c3 = c^3; c4 = c^4; c5 = c^5;
ac = (a-c)/(a+c);                                   % [-]
coef1 = 2*pi*a*b2;                                  % [m^3]
den1 = 8*pi*(1-nu);                                 % [-]

Q = 3/den1;                                         % [-]       - parameter from [1]
R = (1-2*nu)/den1;                                  % [-]       - parameter from [1]
Ia = -coef1*(2/(a*c2) + log(ac)/c3);                % [-]       - parameter from [1]
Iaa = -coef1*(2/(3*a3*c2) + 2/(a*c4) + log(ac)/c5); % [1/m^2]  - parameter from [1]

a11 = 2*R*(Ia-4*pi);                                % [-]        - (A-1) from [2]
a12 = -2*R*(Ia+4*pi);                               % [-]        - (A-2) from [2]
a21 = Q*a2*Iaa + R*Ia - 1;                          % [-]        - (A-3) from [2]
a22 = -Q*a2*Iaa - Ia*(2*R-Q);                       % [-]        - (A-4) from [2]

den2 = 3*lambda+2*mu;                               % [Pa]
num2 = 3*a22-a12;                                   % [-]
den3 = a11*a22-a12*a21;                             % [-]
num3 = a11-3*a21;                                   % [-]

Pdila = P*(2*mu/den2)*(num2-num3)/den3;                     % [units of P]  - (A-5) from [2]
Pstar = P*(1/den2)*(num2*lambda+2*(lambda+mu)*num3)/den3;   % [units of P]  - (A-6) from [2]

a1 = - 2*b2*Pdila;                                  % [m^2*Pa]  - force from [1]
b1 = 3*(b2/c2)*Pdila + 2*(1-2*nu)*Pstar;            % [Pa]      - pressure from [1]
