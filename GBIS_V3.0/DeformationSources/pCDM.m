function [ue,un,uv]=pCDM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,DVx,DVy,...
    DVz,nu)
% pCDM
% calculates the surface displacements associated with a point CDM that is 
% composed of three mutually orthogonal point tensile dislocations in a 
% half-space.
% 
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
% Horizontal coordinates (in EFCS) and depth of the point CDM. The depth
% must be a positive value. X0, Y0 and depth have the same unit as X and Y.
%
% omegaX, omegaY and omegaZ:
% Clockwise rotation angles about X, Y and Z axes, respectively, that 
% specify the orientation of the point CDM in space. The input values must 
% be in degrees.
%
% DVx, DVy and DVz:
% Potencies of the PTDs that before applying the rotations are normal to 
% the X, Y and Z axes, respectively. The potency has the unit of volume 
% (the unit of displacements and CDM semi-axes to the power of 3).
%
% nu:
% Poisson's ratio.
% 
% 
% OUTPUTS
% ue, un and uv:
% Calculated displacement vector components in EFCS. ue, un and uv have the
% same unit as X and Y in inputs.
%
%
% Example-1: Calculate and plot the vertical displacements on a regular 
% grid
%
% [X,Y] = meshgrid(-7:.02:7,-5:.02:5);
% X0 = 0.5; Y0 = -0.25; depth = 2.75; omegaX = 5; omegaY = -8; omegaZ = 30;
% DVx = 0.00144; DVy = 0.00128; DVz = 0.00072; nu = 0.25;
% [ue,un,uv] = pCDM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,DVx,DVy,DVz,nu);
% figure
% surf(X,Y,reshape(uv,size(X)),'edgecolor','none')
% view(2)
% axis equal
% axis tight
% set(gcf,'renderer','painters')
% 
% Example-2: Comparing the pCDM with the center of dilatation (CD)
% 
% The surface displacements associated with a center of dilatation with a 
% volume change of dV = 0.001 are equivalent to those of a point CDM with
% potencies DVx = DVy = DVz = 0.0006. These values for the potencies are 
% calculated from equation 10 in the reference journal article (see below).
% Numerical verification:
% 
% X = -10:0.01:10; Y = zeros(size(X)); Z = zeros(size(X)); X0 = 0; Y0 = 0; 
% depth = 3; dV = 0.001; nu = 0.25; 
% 
% omegaX = 20*randn(1); omegaY = 10*randn(1); omegaZ = 30*randn(1);
% DVx = 0.0006; DVy = 0.0006; DVz = 0.0006;
% 
% [ue1,~,uv1,DV1] = CDdisp(X,Y,Z,X0,Y0,depth,dV,nu);
% [ue2,~,uv2] = pCDM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,DVx,DVy,DVz,nu);
% figure
% plot(X,ue1/max(uv1),'r','LineWidth',2)
% hold on
% plot(X(1:30:end),ue2(1:30:end)/max(uv2),'k.','MarkerSize',15)
% plot(X,uv1/max(uv1),'r','LineWidth',2)
% plot(X(1:30:end),uv2(1:30:end)/max(uv2),'k.','MarkerSize',15)
% xlabel('X')
% ylabel('Normalized displacements')
% text(7,0.2,'Radial')
% text(2,0.75,'Vertical')
% legend('CD','pCDM')
% 
% Example-3: Comparing the pCDM with the CDM (in the far field)
% 
% The surface displacements associated with a CDM that is in the far field 
% are very similar to those of a point CDM with the same coordinates and 
% rotation angles, and with potencies of DVx = 4*ay*az*opening,
% DVy = 4*ax*az*opening and DVz = 4*ax*ay*opening. For more details check 
% the reference journal article (see below).
% Numerical verification:
% 
% X = -10:0.01:10; Y = zeros(size(X)); X0 = 0; Y0 = 0; depth = 4; 
% omegaX = 0; omegaY = 0; omegaZ = 0; ax = 1.25; ay = 1.0;
% az = 0.35; opening = 0.005; nu = 0.25; 
% DVx = 4*ay*az*opening; DVy = 4*ax*az*opening; DVz = 4*ax*ay*opening;
% 
% [ue1,~,uv1] = CDM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,ax,ay,az,...
% opening,nu);
% [ue2,~,uv2] = pCDM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,DVx,DVy,DVz,nu);
% figure
% plot(X,ue1/max(uv1),'r','LineWidth',2)
% hold on
% plot(X(1:30:end),ue2(1:30:end)/max(uv2),'k.','MarkerSize',15)
% plot(X,uv1/max(uv1),'r','LineWidth',2)
% plot(X(1:30:end),uv2(1:30:end)/max(uv2),'k.','MarkerSize',15)
% xlabel('X')
% ylabel('Normalized displacements')
% text(6,0.2,'Horizontal')
% text(2,0.75,'Vertical')
% legend('CDM','point CDM')

% Reference journal article:
% Nikkhoo, M., Walter, T. R., Lundgren, P. R., Prats-Iraola, P. (2016):
% Compound dislocation models (CDMs) for volcano deformation analyses.
% Submitted to Geophysical Journal International

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
% Created: 2015.5.22
% Last modified: 2016.10.18
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

% The potencies must have the same sign
DVsign = [sign(DVx) sign(DVy) sign(DVz)];
if any(DVsign>0) && any(DVsign<0)
    error('Input error: DVx, DVy and DVz must have the same sign!')
end

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
