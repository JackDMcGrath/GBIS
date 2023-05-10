function [m_out]=fault_patches(m,ns,nd,plot_flag)
%FAULT_PATCHES divide a fault into patches
%   [M_OUT]=FAULT_PATCHES(M,NS,ND)
%
%   M is a 7x1 vector based on the disloc parameters
%	   LENGTH : length in the STRIKE direction (LENGTH > 0)
%	   WIDTH  : width in the DIP direction (WIDTH > 0)
%	   DEPTH  : depth to top (DIP<0) or bottom edge (DIP>0) (DEPTH > 0)
%	   DIP    : angle between the fault and a horizontal plane (0 to 90 deg)
%	   STRIKE : dislocation strike 
%      X0     : X-coordinate of centre of edge at depth DEPTH
%      Y0     : Y-coordinate of centre of edge at depth DEPTH
%   NS is number of patches along strike
%   ND is number of patches down dip
%   PLOT_FLAG is 1 (default) for plot, 0 otherwise 
%
%   Andy Hooper, 2022

if nargin<3
    help fault_patches
    error('not enough parameters')
end

if nargin<4
    plot_flag=1;
end

    
dip=m(4)*pi/180;
strike=m(5)*pi/180;

N=ns*nd;
pixs=m(1)/ns; % pixels size along strike
pixd=m(2)/nd; % pixels size along dip

[AS,UD]=meshgrid([0:pixs:(ns-1)*pixs]-(ns-1)*pixs/2,[0:pixd:(nd-1)*pixd]);
AS=AS(:)'; % along strike
UD=UD(:)'; % up dip

m_out=zeros(10,N);
m_out(1,:)=repmat(pixs,1,N);
m_out(2,:)=repmat(pixd,1,N);
m_out(3,:)=repmat(m(3),1,N)-UD*sin(dip);
m_out(4,:)=repmat(m(4),1,N);
m_out(5,:)=repmat(m(5),1,N);
m_out(6,:)=AS*sin(strike)-UD*cos(strike)*cos(dip)+m(6);
m_out(7,:)=AS*cos(strike)+UD*sin(strike)*cos(dip)+m(7);

if plot_flag
    figure
    drawmodel(m_out)
    axis equal;
end



