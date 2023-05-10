function P = P(h,x)
% INPUT
%
% h     dimensionless source depth
% x     dummy variable
%
% OUTPUT
% P1-P4 are expressions from Appendix A of Fialko et al (2001). Used in the 
% definition of the functions T1 - T4 (see T.m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fialko, Y, Khazan, Y and M. Simons (2001). Deformation due to a 
% pressurized horizontal circular crack in an elastic half-space, with 
% applications to volcano geodesy. Geophys. J. Int., 146, 181–190
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


P(1,:) = (12*h^2-x.^2)./(4*h^2+x.^2).^3;                                    % P1(x), pg 189
P(2,:) = log(4*h^2+x.^2) + (8*h^4+2*x.^2*h^2-x.^4)./(4*h^2+x.^2).^2;        % P2(x), pg 189
P(3,:) = 2*(8*h^4-2*x.^2*h^2+x.^4)./(4*h^2+x.^2).^3;                        % P3(x), pg 189
P(4,:) = (4*h^2-x.^2)./(4*h^2+x.^2).^2;                                     % P4(x), pg 189


