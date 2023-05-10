function [V]=arctan(dist,arcpar,offsetVel)
% Now using variable offset based of AUS fixed (taking AUS rate as being
% ~ 97.5 % of plate motion on that side of the fault (8 pi * LD)
% Arcpar = [Velocity, Locking Depth, Rigidity, Horizontal Offset]

if nargin < 3
    offsetVel = 1;
end

    ix=find(dist>=arcpar(4));
    iy=find(dist<arcpar(4));
    
    %% Vertical Offset
    total_vel=-arcpar(1)*(1-arcpar(3));
    ausdist=-8*pi*arcpar(2);
    if offsetVel == 1
        ausvel=((2*(1-arcpar(3))*arcpar(1))/pi).*atan(ausdist/arcpar(2)); % No need to consider a horizontal offset for this
    else
        ausvel = 0;
    end
%     fprintf(' Assym: %.1f\n Total Vel: %.1f mm/yr\n Aus Dist: %.1f km\n Aus Vel: %.1f mm/yr\n\n',arcpar(3),total_vel*1e3,ausdist*1e-3,ausvel*1e3)
    V(ix)=((2*arcpar(3)*arcpar(1))/pi).*atan((dist(ix)-arcpar(4))/arcpar(2))-ausvel;
    V(iy)=(2*(1-arcpar(3))*arcpar(1)/pi).*atan((dist(iy)-arcpar(4))/arcpar(2))-ausvel;
end