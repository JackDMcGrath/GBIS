function [T1 T2 T3 T4] = T(h,t,r)
% INPUT
%
% h     dimensionless source depth
% r     dummy variable (used to integrate along the sill dimenionless radius) 
% t     dummy variable
%
% OUTPUT
% T1-T4 are expressions from Appendix A of Fialko et al (2001), equations 
% (A2) - (A5). T1-T4 are used in the definition of the functions phi and 
% psi (see phi_psi.m)
%
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


M = length(t); N = length(r); 
T1 = zeros(M,N); T2 = T1; T3 = T1;                                          % pre-allocate variables

for i=1:M
    Pm = P(h,t(i)-r);                                                       % define functions P1-P4
    Pp = P(h,t(i)+r);                                                       % define functions P1-P4
    T1(i,:) = 4*h^3*(Pm(1,:)-Pp(1,:));                                      % equation (A2)
    T2(i,:) = (h./(t(i).*r)).*(Pm(2,:)-Pp(2,:)) +h*(Pm(3,:)+Pp(3,:));       % equation (A3)
    T3(i,:) = (h^2./r).*(Pm(4,:)-Pp(4,:)...
                              -2*r.*((t(i)-r).*Pm(1,:)+(t(i)+r).*Pp(1,:))); % equation (A4)
end;
    T4 = T3';                                                               % equation (A5)
