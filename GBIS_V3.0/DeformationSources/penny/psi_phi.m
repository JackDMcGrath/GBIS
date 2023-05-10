function [phi psi t w] = psi_phi(h)
% function [phi psi] = psi_phi(h,t)
% compute the function phi(t) and psi(t) using the Nystrom routine with the
% N-point Gauss-Legendre rule (Numerical Recipes in Fortran 77 (1992), 18.1
% Fredholm Equations of the Second Kind, p. 782.)
%
% INPUT
% h     dimensionless source depth
% t     dummy variable vector (value between 0 and 1; length(t) = M)
%
% OUTPUT
% phi and psi are expressions from Appendix A of Fialko et al (2001), 
% equation (A1). phi and psi are used in the definition of the functions 
% PHI and PSI. See also A_B.m and Fialko et al (2001), eq. (26), p. 184.
% t and w are the abscissas and weights from the Gauss-Legendre quadrature
% rule
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fialko, Y, Khazan, Y and M. Simons (2001). Deformation due to a 
%  pressurized horizontal circular crack in an elastic half-space, with 
%  applications to volcano geodesy. Geophys. J. Int., 146, 181–190
% Press W. et al. (1992). Numerical Recipes in Fortran 77. 2nd Ed, 
%  Cambridge University Press. Available on-line at www.nrbook.com/fortran/ 
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


[t w] = gauleg(0,1,41);                                                     % abscissas and weights from the Gauss-Legendre quadrature rule

% Solution at the quadrature points r(j) **********************************
g = -2*t/pi;                                                                % see eq. (A1) in Fialko et al. (2001)
d = [g zeros(size(g))]';                                                     

[T1 T2 T3 T4] = T(h,t,t);                                                   % Fredholm's integration kernels (Falko et al., 2001, eq. 27)
T1tilde = zeros(size(T1)); T2tilde = zeros(size(T1));                       % pre-allocate variable
T3tilde = zeros(size(T1)); T4tilde = zeros(size(T1));

N = length(t); 
for j= 1:N                                                                  % multiply the integration kernels by the quadrature weights
    T1tilde(:,j) = w(j)*T1(:,j);                                            
    T2tilde(:,j) = w(j)*T2(:,j);
    T3tilde(:,j) = w(j)*T3(:,j);
    T4tilde(:,j) = w(j)*T4(:,j);
end;

Ktilde = [T1tilde, T3tilde; T4tilde, T2tilde];                              % see eq (18.1.5) in Press et al (1992)
y = (eye(2*N,2*N)-(2/pi)*Ktilde)\d;                                         % solution of linear system, (18.1.6) in Press et al (1992)
phi = y(1:N);                                                               % phi at the quadrature points t(j)
psi = y(N+1:2*N);                                                           % psi at the quadrature points t(j)
% *************************************************************************