function [Xnew,Ynew]=point_to_lineXY(X,Y,line)
% Brief function to take a point XY, and move it to the nearest point of a
% line (given as a function of co-ords in 2*N vector)

d=sqrt(sum((line-[X;Y]).^2)); % Find distance to each point of trace
if numel(find(d==min(d)))==1
    dmin(1)=find(d==min(d));
else % In case there are 2 points the same distance away
    f=find(d==min(d));
    dmin(1)=f(1);
end
d(dmin(1))=d(dmin(1))*1e12; % Make that minimum distance massive to find next smallest
dmin(2)=find(d==min(d)); % Location of ends of the nearest segment

linex=line(1,dmin); % X-coords of nearest segment
liney=line(2,dmin); % Y-coords of nearest segment

coeffs=polyfit(linex,liney,1); % Find y=mx+c of line segment
m=coeffs(1);
c1=coeffs(2);


% % Old Version - tried to calculate fault normal to the point, to plot the
% % along the fault segment. Currently glitchy (causing GBIS to get stuck on
% % one spot)
% c2=Y+X/m; % Rearrange y=mx+c to find c for line normal to the segment, passing through XY
% 
% Xnew=(m*(c2-c1))/(m^2+1); %X-coordinate of intersection of interestion using segment mx+c
% Ynew=m*Xnew+c1; % Y-coordinate of interestion using segment mx+c
% 
% Y_check=-(Xnew/m)+c2; % Y-coordinate of interestion using normal mx+c
% 
% if round(Ynew)~=round(Y_check) % Check Y's are the same, round to account for numerical instability that sometimes occurs with polyfit
% fprintf('Error: Inconsistent Y-Coordinate found during check\n')
% clear Xnew Ynew
% end


% New, less idea version. Keep X-coordinate the same, just change Y to make
% it plot along a projection of the nearest fault segment

Xnew = X;
Ynew = m * Xnew + c1;