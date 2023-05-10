function [d,s] = point_to_line(pt, V1, V2)
% Function to find the shortest distance to a line from a point, using the 
% edges of the line.
% Taken from mathworks website.
% pt should be nx3
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
% v1 = [0,0,0];
% v2 = [3,0,0];
% pt = [0,5,0];
% distance = point_to_line(pt,v1,v2)
if size(pt,2)==2
    pt=[pt,zeros(size(pt,1),1)];
end

if size(V1,2)==2
    V1=[V1,0];
end

if size(V2,2)==2
    V2=[V2,0];
end

v1 = repmat(V1,size(pt,1),1);
v2 = repmat(V2,size(pt,1),1);
a = v1 - v2;
b = pt - v2;
d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));

% To work out which side of a line it is on
% https://math.stackexchange.com/questions/274712/calculate-on-which-side-of-a-straight-line-is-a-given-point-located

s=(pt(:,1)-V1(1)).*(V2(2)-V1(2))-(pt(:,2)-V1(2)).*(V2(1)-V1(1));
s=s./abs(s); % sets one side to -1 and the other to 1
ix=isnan(s);
s(ix)=1; % sets NaNs (zero distance) to one
d=s.*d;

end