function [dist] = dist2fault_trace(fault, obs)
% Function to take the UTM position of a pixel, and find the distance to a
% UTM fault trace


    xlimit=[min(obs(1,:)) max(obs(1,:))];
    ylimit=[min(obs(2,:)) max(obs(2,:))];
    
    xbox=xlimit([1 1 2 2 1]);
    ybox=ylimit([1 2 2 1 1]);
    
    % Find where fault intersects with the bounding box
    [xi,yi]=polyxpoly(fault(:,1),fault(:,2),xbox,ybox);
    
    % Find fault segment that crosses the data CAUSING NANS
%     fault1=find(fault(:,1)>min(xi));
%     fault1=fault1(1);
%     fault2=find(fault(:,1)<max(xi));
%     fault2=fault2(end);
%     fault_section=fault([fault1 fault2],:);
%     [XX,YY]=meshgrid(obs(1,:),obs(2,:));
%     [dist] = point_to_line(obs', fault_section(1,:), fault_section(2,:));
    
    % Just work out distance based off straight line between where the
    % fault intersects the bounding box (assumes first two interesctions in
    % the event of there being more than 2)
%     [dist] = point_to_line(obs', xi(1:2)', yi(1:2)');
    [dist] = point_to_line(obs', [xi(1),yi(1)], [xi(2),yi(2)]);