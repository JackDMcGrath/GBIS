function [x w] = gauleg(x1,x2,N)
% function [x w] = gauleg(x1,x2,N)
% Given the upper and lower limits of integration x1 and x2, and given N,
% this routine returns arrays x(1:n) and w(1:n) of length N, containing the
% abscissas x and weights w of the Gaussian- legendre n-point quadrature
% formula. 
% *************************************************************************
% References.
% Loosely based on SUBROUTINE gauleg (Numerical Recipes in Fortran 77, 4.5) 
% *************************************************************************
% Note
% tested against the "High-precision Abscissas and Weights for Gaussian
% Quadrature of high order." table from http://www.holoborodko.com/pavel/ 
% click on Numerical Methods then Numerical Integration. The Table is at
% the end of the page
% ****** for N = 11
%  Abscissas x                  Weights w
%  0                            0.2729250867779006307144835
% ±0.2695431559523449723315320 	0.2628045445102466621806889
% ±0.5190961292068118159257257 	0.2331937645919904799185237
% ±0.7301520055740493240934163 	0.1862902109277342514260976
% ±0.8870625997680952990751578 	0.1255803694649046246346943
% ±0.9782286581460569928039380 	0.0556685671161736664827537
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


z = zeros(1,N);         % pre-allocate variable
xm = 0.5*(x2+x1);       % mid-point
xl = 0.5*(x2-x1);       % half-interval

for n=1:N               % loop over the desidered roots
    z(n) = cos(pi*(n-0.25)/(N+0.5));          % approximation of the ith root of the Legendre's polynomials
    z1 = 100*z(n);
    while abs(z1-z(n)) > eps                   % Newton's method
        [pN, dpN] = legpol(z(n),N+1);          % compute the Legendre's polynomial and its derivative
        z1 = z(n);
        z(n) = z1 - pN/dpN;
    end;
end;    
    [~, dpN] = legpol(z,N+1);               % compute the derivative of the Legendre's polynomial
    x(1:N) = xm - xl*z;                     % scale the root to the desidered interval
    w(1:N) = 2*xl./((1-z.^2).*dpN.^2);      % compute the weights

    