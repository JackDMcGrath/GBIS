function [ue,un,uv,dV,DV]=pECM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,...
    ax,ay,az,p,mu,lambda)
% pECM
% calculates the surface displacements caused by a uniformly-pressurized
% point ellipsoidal cavity in a uniform elastic half-space.
% 
% pECM: point Ellipsoidal Cavity Model
% pCDM: point Compound Dislocation Model
% PTD: Point Tensile Dislocation
% EFCS: Earth-Fixed Coordinate System
%
% INPUTS
% X and Y:
% Horizontal coordinates of calculation points in EFCS (East, North, Up).
% X and Y must have the same size.
%
% X0, Y0 and depth:
% Horizontal coordinates (in EFCS) and depth of the point ECM. The depth
% must be a positive value. X0, Y0 and depth have the same unit as X, Y and
% Z.
%
% omegaX, omegaY and omegaZ:
% Clockwise rotation angles about X, Y and Z axes, respectively, that 
% specify the orientation of the point ECM in space. The input values must 
% be in degrees.
%
% ax, ay and az:
% Semi-axes of the pECM along the X, Y and Z axes, respectively, before
% applying the rotations. ax, ay and az have the same unit as X and Y.
%
% p:
% Pressure on the cavity walls. p has the same unit as the Lamé constants. 
% 
% mu and lambda:
% Lamé constants.
% 
% 
% OUTPUTS
% ue, un and uv:
% Calculated displacement vector components in EFCS. ue, un and uv have the
% same unit as X, Y and Z in inputs.
%
% dV and DV:
% Volume change and potency of the point ECM. The units of the volume 
% change and potency are the same (the unit of displacements and point ECM 
% semi-axes to the power of 3).
% 
% 
% Example: Calculate and plot the vertical displacements on a regular grid
%
% [X,Y] = meshgrid(-7:.02:7,-5:.02:5);
% X0 = 0.5; Y0 = -0.25; depth = 2.75; omegaX = 5; omegaY = -8; omegaZ = 30;
% ax = 1; ay = 0.75; az = 0.25; p = 1e6; mu = 0.2e10; lambda = 0.2e10;
% [~,~,uv,dV,DV] = pECM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,...
%     ax,ay,az,p,mu,lambda);
% figure
% surf(X,Y,reshape(uv,size(X)),'edgecolor','none')
% view(2)
% axis equal
% axis tight
% set(gcf,'renderer','painters')

% Reference journal articles:
% 1)
% Nikkhoo, M., Walter, T. R., Lundgren, P. R., Prats-Iraola, P. (2017):
% Compound dislocation models (CDMs) for volcano deformation analyses.
% Submitted to Geophysical Journal International, doi: 10.1093/gji/ggw427
% 
% 2)
% Eshelby, J. D. (1957):
% The determination of the elastic field of an ellipsoidal inclusion, and
% related problems.
% Proceedings of the royal society of London. Series A. Mathematical and
% physical sciences. 241 (1226), 376-396. doi: 10.1098/rspa.1957.0133
%
% 3)
% Carlson, B. C. (1995):
% Numerical computation of real or complex elliptic integrals.
% Numer. Algor., 10(1), 13–26. doi: 10.1007/BF02198293
% 
% 4)
% Okada, Y. (1992):
% Internal deformation due to shear and tensile faults in a half-space, 
% Bull. seism. Soc. Am., 82(2), 1018–1040.

% Copyright (c) 2016 Mehdi Nikkhoo
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

% I appreciate any comments or bug reports.

% Mehdi Nikkhoo
% Created: 2015.7.24
% Last modified: 2022.10.24
%
% Section 2.1, Physics of Earthquakes and Volcanoes
% Department 2, Geophysics
% Helmholtz Centre Potsdam
% German Research Centre for Geosciences (GFZ)
%
% email:
% mehdi.nikkhoo@gfz-potsdam.de
% mehdi.nikkhoo@gmail.com
%
% website:
% http://www.volcanodeformation.com

nu = lambda/(lambda+mu)/2; % Poisson's ratio
K = lambda+2*mu/3; % Bulk modulus

% Set a threshold to avoid geometries with zero volume
r0 = 1e-12; % Threshold for stability in the shape tensor calculations
ax(ax<r0) = r0;
ay(ay<r0) = r0;
az(az<r0) = r0;
[ai,Ind] = sort([ax ay az],'descend');

S = ShapeTensorECM(ai(1),ai(2),ai(3),nu);
Sm = [S(1)-1 S(2) S(3);S(4) S(5)-1 S(6);S(7) S(8) S(9)-1]; % Shape tensor

eT = -inv(Sm)*p*ones(3,1)/3/K; % Transformation strain
% For uniformly-pressurized ellipsoids the eT elements have the same sign!
eT(sign(eT)~=sign(p)) = 0;
V = 4/3*pi*ax*ay*az; % Ellipsoid volume

% calculate actual volume change
dV = (sum(eT)-p/K)*V;
% calculate potency
DV = dV+p*V/K; % Also DV = sum(M)/3/K could be used!

tmp = zeros(3,1);
tmp(Ind) = V*eT;
DVx = tmp(1);
DVy = tmp(2);
DVz = tmp(3);

% The point CDM calculation
X = X(:);
Y = Y(:);

Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

Vstrike1 = [-R(2,1),R(1,1),0];
Vstrike1 = Vstrike1/norm(Vstrike1);
strike1 = atan2(Vstrike1(1),Vstrike1(2))*180/pi;
if isnan(strike1)
    strike1 = 0;
end
dip1 = acosd(R(3,1));

Vstrike2 = [-R(2,2),R(1,2),0];
Vstrike2 = Vstrike2/norm(Vstrike2);
strike2 = atan2(Vstrike2(1),Vstrike2(2))*180/pi;
if isnan(strike2)
    strike2 = 0;
end
dip2 = acosd(R(3,2));

Vstrike3 = [-R(2,3),R(1,3),0];
Vstrike3 = Vstrike3/norm(Vstrike3);
strike3 = atan2(Vstrike3(1),Vstrike3(2))*180/pi;
if isnan(strike3)
    strike3 = 0;
end
dip3 = acosd(R(3,3));

% Calculate contribution of the first PTD
if DVx~=0
    [ue1,un1,uv1] = PTDdispSurf(X,Y,X0,Y0,depth,strike1,dip1,DVx,nu);
else
    ue1 = zeros(size(X));
    un1 = zeros(size(X));
    uv1 = zeros(size(X));
end

% Calculate contribution of the second PTD
if DVy~=0
    [ue2,un2,uv2] = PTDdispSurf(X,Y,X0,Y0,depth,strike2,dip2,DVy,nu);
else
    ue2 = zeros(size(X));
    un2 = zeros(size(X));
    uv2 = zeros(size(X));
end

% Calculate contribution of the third PTD
if DVz~=0
    [ue3,un3,uv3] = PTDdispSurf(X,Y,X0,Y0,depth,strike3,dip3,DVz,nu);
else
    ue3 = zeros(size(X));
    un3 = zeros(size(X));
    uv3 = zeros(size(X));
end

ue = ue1+ue2+ue3;
un = un1+un2+un3;
uv = uv1+uv2+uv3;


function [S]=ShapeTensorECM(a1,a2,a3,nu)
% ShapeTensorECM
% calculates the Eshelby (1957) shape tensor components.

if all([a1 a2 a3]==0)
    S = zeros(9,1);
    return
end

% Calculate Ik and Iij terms for triaxial, oblate and prolate ellipsoids
if a1>a2 && a2>a3 && a3>0
    % General case: triaxial ellipsoid
    sin_theta = sqrt(1-a3^2/a1^2);
    k = sqrt((a1^2-a2^2)/(a1^2-a3^2));
    
    % % Calculate Legendre's incomplete elliptic integrals of the first and
    % % second kind using MATLAB Symbolic Math Toolbox
    % F =  mfun('EllipticF',sin_theta,k);
    % E =  mfun('EllipticE',sin_theta,k);
    
    % Calculate Legendre's incomplete elliptic integrals of the first and
    % second kind using Carlson (1995) method (see Numerical computation of
    % real or complex elliptic integrals. Carlson, B.C. Numerical 
    % Algorithms (1995) 10: 13. doi:10.1007/BF02198293)
    tol = 1e-16;
    c = 1/sin_theta^2;
    F = RF(c-1,c-k^2,c,tol);
    E = F-k^2/3*RD(c-1,c-k^2,c,tol);
    
    I1 = 4*pi*a1*a2*a3/(a1^2-a2^2)/sqrt(a1^2-a3^2)*(F-E);
    I3 = 4*pi*a1*a2*a3/(a2^2-a3^2)/sqrt(a1^2-a3^2)*...
        (a2*sqrt(a1^2-a3^2)/a1/a3-E);
    I2 = 4*pi-I1-I3;
    
    I12 = (I2-I1)/(a1^2-a2^2);
    I13 = (I3-I1)/(a1^2-a3^2);
    I11 = (4*pi/a1^2-I12-I13)/3;
    
    I23 = (I3-I2)/(a2^2-a3^2);
    I21 = I12;
    I22 = (4*pi/a2^2-I23-I21)/3;
    
    I31 = I13;
    I32 = I23;
    I33 = (4*pi/a3^2-I31-I32)/3;
    
elseif a1==a2 && a2>a3 && a3>0
    % Special case-1: Oblate ellipsoid
    I1 = 2*pi*a1*a2*a3/(a1^2-a3^2)^1.5*(acos(a3/a1)-a3/a1*...
        sqrt(1-a3^2/a1^2));
    I2 = I1;
    I3 = 4*pi-2*I1;
    
    I13 = (I3-I1)/(a1^2-a3^2);
    I11 = pi/a1^2-I13/4;
    I12 = I11;
    
    I23 = I13;
    I22 = pi/a2^2-I23/4;
    I21 = I12;
    
    I31 = I13;
    I32 = I23;
    I33 = (4*pi/a3^2-2*I31)/3;
    
elseif a1>a2 && a2==a3 && a3>0
    % Special case-2: Prolate ellipsoid
    I2 = 2*pi*a1*a2*a3/(a1^2-a3^2)^1.5*(a1/a3*sqrt(a1^2/a3^2-1)-...
        acosh(a1/a3));
    I3 = I2;
    I1 = 4*pi-2*I2;
    
    I12 = (I2-I1)/(a1^2-a2^2);
    I13 = I12;
    I11 = (4*pi/a1^2-2*I12)/3;
    
    I21 = I12;
    I22 = pi/a2^2-I21/4;
    I23 = I22;
    
    I32 = I23;
    I31 = I13;
    I33 = (4*pi/a3^2-I31-I32)/3;
end

% Calculate the shape-tensor components
if a1==a2 && a2==a3
    % Special case-3: Sphere
    S1111 = (7-5*nu)/15/(1-nu);
    S1122 = (5*nu-1)/15/(1-nu);
    S1133 = (5*nu-1)/15/(1-nu);
    S2211 = (5*nu-1)/15/(1-nu);
    S2222 = (7-5*nu)/15/(1-nu);
    S2233 = (5*nu-1)/15/(1-nu);
    S3311 = (5*nu-1)/15/(1-nu);
    S3322 = (5*nu-1)/15/(1-nu);
    S3333 = (7-5*nu)/15/(1-nu);
else
    % General triaxial, oblate and prolate ellipsoids
    S1111 = 3/8/pi/(1-nu)*a1^2*I11+(1-2*nu)/8/pi/(1-nu)*I1;
    S1122 = 1/8/pi/(1-nu)*a2^2*I12-(1-2*nu)/8/pi/(1-nu)*I1;
    S1133 = 1/8/pi/(1-nu)*a3^2*I13-(1-2*nu)/8/pi/(1-nu)*I1;
    S2211 = 1/8/pi/(1-nu)*a1^2*I21-(1-2*nu)/8/pi/(1-nu)*I2;
    S2222 = 3/8/pi/(1-nu)*a2^2*I22+(1-2*nu)/8/pi/(1-nu)*I2;
    S2233 = 1/8/pi/(1-nu)*a3^2*I23-(1-2*nu)/8/pi/(1-nu)*I2;
    S3311 = 1/8/pi/(1-nu)*a1^2*I31-(1-2*nu)/8/pi/(1-nu)*I3;
    S3322 = 1/8/pi/(1-nu)*a2^2*I32-(1-2*nu)/8/pi/(1-nu)*I3;
    S3333 = 3/8/pi/(1-nu)*a3^2*I33+(1-2*nu)/8/pi/(1-nu)*I3;
end

S = [S1111 S1122 S1133 S2211 S2222 S2233 S3311 S3322 S3333]';

function [rf]=RF(x,y,z,r)
% RF
% calculates the RF term in the Carlson (1995) method for calculating
% elliptic integrals

% r = 1e-16

if any([x,y,z]<0)
    error('x, y and z values must be positive!')
elseif nnz([x,y,z])<2
    error('At most one of the x, y and z values can be zero!')
end

xm = x;
ym = y;
zm = z;
A0 = (x+y+z)/3;
Q = max([abs(A0-x),abs(A0-y),abs(A0-z)])/(3*r)^(1/6);
n = 0;
Am = A0;
while abs(Am)<=Q/(4^n)
    lambdam = sqrt(xm*ym)+sqrt(xm*zm)+sqrt(ym*zm);
    Am = (Am+lambdam)/4;
    xm = (xm+lambdam)/4;
    ym = (ym+lambdam)/4;
    zm = (zm+lambdam)/4;
    n = n+1;
end
X = (A0-x)/4^n/Am;
Y = (A0-y)/4^n/Am;
Z = -X-Y;
E2 = X*Y-Z^2;
E3 = X*Y*Z;
rf = (1-E2/10+E3/14+E2^2/24-3*E2*E3/44)/sqrt(Am);

function [rd]=RD(x,y,z,r)
% RD
% calculates the RD term in the Carlson (1995) method for calculating
% elliptic integrals

% r = 1e-16

if z==0
    error('z value must be nonzero!')
elseif all([x,y]==0)
    error('At most one of the x and y values can be zero!')
end

xm = x;
ym = y;
zm = z;
A0 = (x+y+3*z)/5;
Q = max([abs(A0-x),abs(A0-y),abs(A0-z)])/(r/4)^(1/6);
n = 0;
Am = A0;
S = 0;
while abs(Am)<=Q/(4^n)
    lambdam = sqrt(xm*ym)+sqrt(xm*zm)+sqrt(ym*zm);
    S = S+(1/4^n)/sqrt(zm)/(zm+lambdam);
    Am = (Am+lambdam)/4;
    xm = (xm+lambdam)/4;
    ym = (ym+lambdam)/4;
    zm = (zm+lambdam)/4;
    n = n+1;
end

X = (A0-x)/4^n/Am;
Y = (A0-y)/4^n/Am;
Z = -(X+Y)/3;
E2 = X*Y-6*Z^2;
E3 = (3*X*Y-8*Z^2)*Z;
E4 = 3*(X*Y-Z^2)*Z^2;
E5 = X*Y*Z^3;
rd = (1-3*E2/14+E3/6+9*E2^2/88-3*E4/22-9*E2*E3/52+3*E5/26)/4^n/Am^1.5+3*S;

function [ue,un,uv]=PTDdispSurf(X,Y,X0,Y0,depth,strike,dip,DV,nu)
% PTDdispSurf calculates surface displacements associated with a tensile 
% point dislocation (PTD) in an elastic half-space (Okada, 1985).

x = X(:)-X0;
y = Y(:)-Y0;

beta = strike-90;
Rz = [cosd(beta) -sind(beta);sind(beta) cosd(beta)];
r_beta = Rz*[x y]';
x = r_beta(1,:)';
y = r_beta(2,:)';

r = (x.^2+y.^2+depth.^2).^0.5;
d = depth;
q = y*sind(dip)-d*cosd(dip);

I1 = (1-2*nu)*y.*(1./r./(r+d).^2-x.^2.*(3*r+d)./r.^3./(r+d).^3);
I2 = (1-2*nu)*x.*(1./r./(r+d).^2-y.^2.*(3*r+d)./r.^3./(r+d).^3);
I3 = (1-2*nu)*x./r.^3-I2;
I5 = (1-2*nu)*(1./r./(r+d)-x.^2.*(2*r+d)./r.^3./(r+d).^2);

% Note: For a PTD M0 = DV*mu!
ue = DV/2/pi*(3*x.*q.^2./r.^5-I3*sind(dip)^2);
un = DV/2/pi*(3*y.*q.^2./r.^5-I1*sind(dip)^2);
uv = DV/2/pi*(3*d.*q.^2./r.^5-I5*sind(dip)^2);

r_beta = Rz'*[ue un]';
ue = r_beta(1,:)';
un = r_beta(2,:)';
