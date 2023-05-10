function [X_Fault,Y_Fault]=faulttrace2faultcenter(m)
% Function to take a position on the fault trace, and to project that to
% the center of the top edge of the corresponding fault

% Coordinates of fault trace
X_Trace = m(6);
Y_Trace = m(7);
X0 = 0;
Y0 = 0;

Z1 = m(3);
dip1 = m(4);
Strike = m(5);

X1 = X0-Z1/tand(dip1);
Y1 = Y0;

% Create vector of Xs and Ys
xs = [X0 X1];
ys = [Y0 Y1];

% Rotate coordinates
xRot = xs*cosd(Strike) - ys*sind(Strike);
yRot = xs*sind(Strike) + ys*cosd(Strike);


X_Fault = xRot(2)+X_Trace;
Y_Fault = -yRot(2)+Y_Trace;

end