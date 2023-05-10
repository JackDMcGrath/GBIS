function [pN dpN] = legpol(x,N)
% legendre polynomial, Bonnet’s recursion formula
%
% Reference 
% http://en.wikipedia.org/wiki/Legendre_polynomials
%
% Note
% index n in wiki goes from 0 to N, in MATLAB j goes from 1 to N+1. To
% obtain the right coefficients we introduced j = n+1 -> n = j-1
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


    P(1,:) = ones(size(x)); dP(1,:) = zeros(size(x));
    P(2,:) = x; 
    for j=2:N                                                               % loop up the recursion relation 
        P(j+1,:) = ((2*j-1)*x.*P(j,:) - (j-1)*P(j-1,:))/j;
         dP(j,:) = (j-1)*(x.*P(j,:) - P(j-1,:))./(x.^2-1);
    end;
    
    pN =  P(N,:);
   dpN = dP(N,:); 