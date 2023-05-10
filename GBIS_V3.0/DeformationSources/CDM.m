function [ue,un,uv,DV]=CDM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,ax,ay,...
    az,opening,nu)
% CDM
% calculates the surface displacements and potency associated with a CDM 
% that is composed of three mutually orthogonal rectangular dislocations in
% a half-space.
% 
% CDM: Compound Dislocation Model
% RD: Rectangular Dislocation
% EFCS: Earth-Fixed Coordinate System
% RDCS: Rectangular Dislocation Coordinate System
% ADCS: Angular Dislocation Coordinate System
% (The origin of the RDCS is the RD centroid. The axes of the RDCS are 
% aligned with the strike, dip and normal vectors of the RD, respectively.)
%
% INPUTS
% X and Y:
% Horizontal coordinates of calculation points in EFCS (East, North, Up).
% X and Y must have the same size.
%
% X0, Y0 and depth:
% Horizontal coordinates (in EFCS) and depth of the CDM centroid. The depth
% must be a positive value. X0, Y0 and depth have the same unit as X and Y.
%
% omegaX, omegaY and omegaZ:
% Clockwise rotation angles about X, Y and Z axes, respectively, that 
% specify the orientation of the CDM in space. The input values must be in 
% degrees.
%
% ax, ay and az:
% Semi-axes of the CDM along the X, Y and Z axes, respectively, before
% applying the rotations. ax, ay and az have the same unit as X and Y.
%
% opening:
% The opening (tensile component of the Burgers vector) of the RDs that
% form the CDM. The unit of opening must be the same as the unit of ax, ay
% and az.
%
% nu:
% Poisson's ratio.
% 
% 
% OUTPUTS
% ue, un and uv:
% Calculated displacement vector components in EFCS. ue, un and uv have the
% same unit as opening and the CDM semi-axes in inputs.
%
% DV:
% Potency of the CDM. DV has the unit of volume (the unit of displacements, 
% opening and CDM semi-axes to the power of 3).
% 
% 
% Example: Calculate and plot the vertical displacements on a regular grid.
%
% [X,Y] = meshgrid(-7:.02:7,-5:.02:5);
% X0 = 0.5; Y0 = -0.25; depth = 2.75; omegaX = 5; omegaY = -8; omegaZ = 30;
% ax = 0.4; ay = 0.45; az = 0.8; opening = 1e-3; nu = 0.25;
% [ue,un,uv,DV] = CDM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,ax,ay,az,...
%     opening,nu);
% figure
% surf(X,Y,reshape(uv,size(X)),'edgecolor','none')
% view(2)
% axis equal
% axis tight
% set(gcf,'renderer','painters')

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

X = X(:);
Y = Y(:);

% convert the semi-axes to axes
ax = 2*ax;
ay = 2*ay;
az = 2*az;

Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

P0 = [X0 Y0 -depth]'; % The centroid

P1 = P0+ay*R(:,2)/2+az*R(:,3)/2;
P2 = P1-ay*R(:,2);
P3 = P2-az*R(:,3);
P4 = P1-az*R(:,3);

Q1 = P0-ax*R(:,1)/2+az*R(:,3)/2;
Q2 = Q1+ax*R(:,1);
Q3 = Q2-az*R(:,3);
Q4 = Q1-az*R(:,3);

R1 = P0+ax*R(:,1)/2+ay*R(:,2)/2;
R2 = R1-ax*R(:,1);
R3 = R2-ay*R(:,2);
R4 = R1-ay*R(:,2);

VertVec = [P1(3) P2(3) P3(3) P4(3) Q1(3) Q2(3) Q3(3) Q4(3)...
    R1(3) R2(3) R3(3) R4(3)];
% if any(VertVec>0)
%     error('Half-space solution: The CDM must be under the free surface!')
% end

if ax==0 && ay==0 && az==0
    ue = zeros(size(X));
    un = zeros(size(X));
    uv = zeros(size(X));
    
elseif ax==0 && ay~=0 && az~=0
    [ue,un,uv] = RDdispSurf(X,Y,P1,P2,P3,P4,opening,nu);
elseif ax~=0 && ay==0 && az~=0
    [ue,un,uv] = RDdispSurf(X,Y,Q1,Q2,Q3,Q4,opening,nu);
elseif ax~=0 && ay~=0 && az==0
    [ue,un,uv] = RDdispSurf(X,Y,R1,R2,R3,R4,opening,nu);
else
    [ue1,un1,uv1] = RDdispSurf(X,Y,P1,P2,P3,P4,opening,nu);
    [ue2,un2,uv2] = RDdispSurf(X,Y,Q1,Q2,Q3,Q4,opening,nu);
    [ue3,un3,uv3] = RDdispSurf(X,Y,R1,R2,R3,R4,opening,nu);
    ue = ue1+ue2+ue3;
    un = un1+un2+un3;
    uv = uv1+uv2+uv3;
end

% Calculate the CDM potency (aX, aY and aZ were converted to full axes)
DV = (ax*ay+ax*az+ay*az)*opening;


function [ue,un,uv]=RDdispSurf(X,Y,P1,P2,P3,P4,opening,nu)
% RDdispSurf calculates surface displacements associated with a rectangular
% dislocation in an elastic half-space.

bx = opening;

Vnorm = cross(P2-P1,P4-P1);
Vnorm = Vnorm/norm(Vnorm);
bX = bx*Vnorm(1);
bY = bx*Vnorm(2);
bZ = bx*Vnorm(3);

[u1,v1,w1] = AngSetupFSC(X,Y,bX,bY,bZ,P1,P2,nu); % Side P1P2
[u2,v2,w2] = AngSetupFSC(X,Y,bX,bY,bZ,P2,P3,nu); % Side P2P3
[u3,v3,w3] = AngSetupFSC(X,Y,bX,bY,bZ,P3,P4,nu); % Side P3P4
[u4,v4,w4] = AngSetupFSC(X,Y,bX,bY,bZ,P4,P1,nu); % Side P4P1

ue = u1+u2+u3+u4;
un = v1+v2+v3+v4;
uv = w1+w2+w3+w4;

function [X1,X2,X3]=CoordTrans(x1,x2,x3,A)
% CoordTrans transforms the coordinates of the vectors, from
% x1x2x3 coordinate system to X1X2X3 coordinate system. "A" is the
% transformation matrix, whose columns e1,e2 and e3 are the unit base
% vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given
% in X1X2X3. The transpose of A (i.e., A') will transform the coordinates
% from X1X2X3 into x1x2x3.

x1 = x1(:);
x2 = x2(:);
x3 = x3(:);
r = A*[x1';x2';x3'];
X1 = r(1,:)';
X2 = r(2,:)';
X3 = r(3,:)';

function [ue,un,uv]=AngSetupFSC(X,Y,bX,bY,bZ,PA,PB,nu)
% AngSetupSurf calculates the displacements associated with an angular
% dislocation pair on each side of an RD in a half-space.

SideVec = PB-PA;
eZ = [0 0 1]';
beta = acos(-SideVec'*eZ/norm(SideVec));

if abs(beta)<eps || abs(pi-beta)<eps
    ue = zeros(length(X),1);
    un = zeros(length(X),1);
    uv = zeros(length(X),1);
else
    ey1 = [SideVec(1:2);0];
    ey1 = ey1/norm(ey1);
    ey3 = -eZ;
    ey2 = cross(ey3,ey1);
    A = [ey1,ey2,ey3]; % Transformation matrix
    
    % Transform coordinates from EFCS to the first ADCS
    [y1A,y2A,~] = CoordTrans(X-PA(1),Y-PA(2),zeros(length(X),1)-PA(3),A);
    % Transform coordinates from EFCS to the second ADCS
    [y1AB,y2AB,~] = CoordTrans(SideVec(1),SideVec(2),SideVec(3),A);
    y1B = y1A-y1AB;
    y2B = y2A-y2AB;
    
    % Transform slip vector components from EFCS to ADCS
    [b1,b2,b3] = CoordTrans(bX,bY,bZ,A);
    
    % Determine the best artefact-free configuration for the calculation
    % points near the free surface
    I = (beta*y1A)>=0;
    
    % Configuration I
    [v1A(I),v2A(I),v3A(I)] = AngDisDispSurf(y1A(I),y2A(I),...
        -pi+beta,b1,b2,b3,nu,-PA(3));
    [v1B(I),v2B(I),v3B(I)] = AngDisDispSurf(y1B(I),y2B(I),...
        -pi+beta,b1,b2,b3,nu,-PB(3));
    
    % Configuration II
    [v1A(~I),v2A(~I),v3A(~I)] = AngDisDispSurf(y1A(~I),y2A(~I),...
        beta,b1,b2,b3,nu,-PA(3));
    [v1B(~I),v2B(~I),v3B(~I)] = AngDisDispSurf(y1B(~I),y2B(~I),...
        beta,b1,b2,b3,nu,-PB(3));
    
    % Calculate total displacements in ADCS
    v1 = v1B-v1A;
    v2 = v2B-v2A;
    v3 = v3B-v3A;
    
    % Transform total displacements from ADCS to EFCS
    [ue,un,uv] = CoordTrans(v1,v2,v3,A');
end

function [v1,v2,v3] = AngDisDispSurf(y1,y2,beta,b1,b2,b3,nu,a)
% AngDisDispSurf calculates the displacements associated with an angular
% dislocation in a half-space.

sinB = sin(beta);
cosB = cos(beta);
cotB = cot(beta);
z1 = y1*cosB+a*sinB;
z3 = y1*sinB-a*cosB;
r2 = y1.^2+y2.^2+a^2;
r = sqrt(r2);

Fi = 2*atan2(y2,(r+a)*cot(beta/2)-y1); % The Burgers function

v1b1 = b1/2/pi*((1-(1-2*nu)*cotB^2)*Fi+y2./(r+a).*((1-2*nu)*(cotB+y1./...
    2./(r+a))-y1./r)-y2.*(r*sinB-y1)*cosB./r./(r-z3));
v2b1 = b1/2/pi*((1-2*nu)*((.5+cotB^2)*log(r+a)-cotB/sinB*log(r-z3))-1./...
    (r+a).*((1-2*nu)*(y1*cotB-a/2-y2.^2./2./(r+a))+y2.^2./r)+y2.^2*...
    cosB./r./(r-z3));
v3b1 = b1/2/pi*((1-2*nu)*Fi*cotB+y2./(r+a).*(2*nu+a./r)-y2*cosB./...
    (r-z3).*(cosB+a./r));

v1b2 = b2/2/pi*(-(1-2*nu)*((.5-cotB^2)*log(r+a)+cotB^2*cosB*log(r-z3))-...
    1./(r+a).*((1-2*nu)*(y1*cotB+.5*a+y1.^2./2./(r+a))-y1.^2./r)+z1.*(r*...
    sinB-y1)./r./(r-z3));
v2b2 = b2/2/pi*((1+(1-2*nu)*cotB^2)*Fi-y2./(r+a).*((1-2*nu)*(cotB+y1./...
    2./(r+a))-y1./r)-y2.*z1./r./(r-z3));
v3b2 = b2/2/pi*(-(1-2*nu)*cotB*(log(r+a)-cosB*log(r-z3))-y1./(r+a).*(2*...
    nu+a./r)+z1./(r-z3).*(cosB+a./r));

v1b3 = b3/2/pi*(y2.*(r*sinB-y1)*sinB./r./(r-z3));
v2b3 = b3/2/pi*(-y2.^2*sinB./r./(r-z3));
v3b3 = b3/2/pi*(Fi+y2.*(r*cosB+a)*sinB./r./(r-z3));

v1 = v1b1+v1b2+v1b3;
v2 = v2b1+v2b2+v2b3;
v3 = v3b1+v3b2+v3b3;
