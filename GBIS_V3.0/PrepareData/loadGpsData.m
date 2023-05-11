function [gps, obsGps, nObsGps] = loadGpsData(gps, geo,plotfig)

% Function to ingest GPS data from text files
%
% Usage: [gps, obsGps, nObsGps] = loadGpsData(gps, geo)
% Input Parameters:
%       gps: gps structure containing path to data files
%       geo: structure with local coordinates origin and bounding box
%
% Output Parameters:
%       gps: structure with GPS data and further information (e.g., inverse
%       of covariance matrix)
%       obsGps: coordinates of observation points
%       nObsGps: number of observation points
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Markov chain Monte Carlo algorithm incorporating the Metropolis alghoritm
% (e.g., Mosegaard & Tarantola, JGR,(1995).
%
% by Marco Bagnardi and Andrew Hooper (COMET, University of Leeds)
% Email: M.Bagnardi@leeds.ac.uk
% Reference: TBA (Bagnardi and Hooper, in prep.)
%
% The function uses third party software.
% =========================================================================
% Last update: 03/05/2017

if nargin<3
    plotfig=1;
end

global outputDir  % Set global variables

gpsTxt = load(gps.dataPath); % Read text file with GPS data

gps.ll = gpsTxt(:, 1:2); % Read GPS site Longitude and Latitude coordinates 
nGps = size(gpsTxt, 1); % Retrieve number of GPS sites
gps.displacements = gpsTxt(:, [5,3,7])'/1000; % Import GPS displacements in mm and tranforms to m
gps.sigmas = gpsTxt(:, [6,4,8])'/1000; % Import GPS displacement st. dev. in mm and tranforms to m

obsGps = llh2local([gps.ll'; zeros(1,nGps)], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
obsGps = [obsGps; zeros(1,size(obsGps,2))]; % Add zeros to third column of observation matrix
nObsGps = size(obsGps,2); % Determine number of entries in GPS observation matrix
gps.variance = gps.sigmas.^2; % Calculate variance of GPS displacements
gps.invCov = diag(1./reshape(gps.variance(1:3,:),nGps*3,1)); % Generate inverse of covariance matrix for GPS
if isfield(gps,'xsillExp')
    obs = obsGps(1:2,:)';
    [X1,X2] = meshgrid(obs(:,1)); % Create square matrices of Xs
    [Y1,Y2] = meshgrid(obs(:,2)); % Create square matrices of Ys
    H = sqrt((X1-X2).^2 + (Y1 - Y2).^2); % Calculate distance between points
    invCov=zeros(nObsGps*3);
    fprintf('Calculating covariance matrix of GPS\n')
    xcovarianceMatrix = gps.xsillExp * exp(-3*H/gps.xrange) + gps.xnugget*eye(nObsGps); % Calculate covariance matrix for exponential model with nugget
    ycovarianceMatrix = gps.ysillExp * exp(-3*H/gps.yrange) + gps.ynugget*eye(nObsGps); % Calculate covariance matrix for exponential model with nugget
    zcovarianceMatrix = gps.zsillExp * exp(-3*H/gps.zrange) + gps.znugget*eye(nObsGps); % Calculate covariance matrix for exponential model with nugget
    xinvCov = inv(xcovarianceMatrix); % Calculate inverse of covariance matrix
    yinvCov = inv(ycovarianceMatrix); % Calculate inverse of covariance matrix
    zinvCov = inv(zcovarianceMatrix); % Calculate inverse of covariance matrix
    [xi,yi]=meshgrid(1:3:nObsGps);
    invCov(1:3:nObsGps*3,1:3:nObsGps*3)=xinvCov;
    invCov(2:3:nObsGps*3,2:3:nObsGps*3)=yinvCov;
    invCov(3:3:nObsGps*3,3:3:nObsGps*3)=yinvCov;
    invCov(find(gps.invCov))=diag(gps.invCov); % Replace variance with that calculated from the variance of the displacements
    gps.invCov=invCov;
end

% Display Gps vectors
obsGps(:,end+1) = [max(obsGps(1,:))+5000; min(obsGps(2,:))-5000; 0]; % add coordinates of legend
hscalebar = max(round(max(abs(gps.displacements(1:2,:)')),3)); % Determine length of scalebar
vscalebar = round(max(abs(gps.displacements(3,:))),3);
gps.displacements(:,end+1) = [-hscalebar 0 0]; % add displacements for legend

if plotfig==0
    vis='off';
    fprintf('GPS figure not plotted. Change flag to 1\n')
elseif plotfig==1
    vis='on';
end
figure('Visible',vis)
quiver(obsGps(1,:), obsGps(2,:), gps.displacements(1,:), gps.displacements(2,:), 1, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 0.03, 'Marker', 's')
axis equal; 
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from origin (m)')
ylabel('Y distance from origin (m)')
title('GPS horizontal displacements')
xlim([min(obsGps(1,:))-10000 max(obsGps(1,:))+10000]);
ylim([min(obsGps(2,:))-10000 max(obsGps(2,:))+10000]);
text(obsGps(1,end),obsGps(2,end)-1500,[num2str(hscalebar*1000),' mm'], 'HorizontalAlignment','Right')
if strcmpi(vis,'on');drawnow;end
if ~isempty(outputDir)
saveas(gcf,[outputDir,'/Figures/GPS_displacements.png'])
end

figure('Visible',vis)
subplot(1,2,1)
quiver(obsGps(1,:), obsGps(2,:), gps.displacements(1,:), gps.displacements(2,:), 1, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 0.03)
axis equal; 
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from origin (m)')
ylabel('Y distance from origin (m)')
title('GPS horizontal displacements')
xlim([min(obsGps(1,:))-10000 max(obsGps(1,:))+10000]);
ylim([min(obsGps(2,:))-10000 max(obsGps(2,:))+10000]);
text(obsGps(1,end),obsGps(2,end)-2500,[num2str(hscalebar*1000),' mm'], 'HorizontalAlignment','Right')

subplot(1,2,2)
quiver(obsGps(1,:), obsGps(2,:), [zeros(1,size(gps.displacements,2)-1),-vscalebar], gps.displacements(3,:), 1, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 0.03)
axis equal; 
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('GPS vertical displacements')
xlim([min(obsGps(1,:))-10000 max(obsGps(1,:))+10000]);
ylim([min(obsGps(2,:))-10000 max(obsGps(2,:))+10000]);
text(obsGps(1,end),obsGps(2,end)-2500,[num2str(vscalebar*1000),' mm'], 'HorizontalAlignment','Right')


if strcmpi(vis,'on');drawnow;end
if ~isempty(outputDir)
saveas(gcf,[outputDir,'/Figures/GPS_displacements_hor_vert.png'])
end

obsGps(:,end) = []; % remove coordinates of legend
gps.displacements(:,end) = []; % remove displacements for legend




