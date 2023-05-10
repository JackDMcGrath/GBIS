function [] = fitVariogram()

% Function to fit an exponential function to the data isotropic (semi-)variogram
%
% Usage:  fitVariogram(inputFile, wavelength)
%   inputFile:  path and name of *.mat file containing data in the format used by
%               GBIS
%
%   wavelength: wavelength of InSAR data in meters (e.g., 0.056 m for
%               Sentine-1/Envisat/ERS)
%
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
% Last update: 02/05/2017

%%
dbstop if error
% close all
clc
inputFile='/nfs/a285/homes/eejdm/InSAR/combined_inversion/UTM_data/mcmc_profiles/fault_testing/GBIS/Data/3comp_insar_prof2_all.txt';
% Load dataset
disp('Ingesting data to estimate (semi-)variogram ...')
gps = load(inputFile);

% Create colormap
cmapSeismo = colormap_cpt('GMT_seis.cpt', 100);    % GMT 'Seismo' colormap for wrapped data

% Find a local reference point
refPoint = [min(gps(:,1)), min(gps(:,2))]; % Determine local reference point

% Convert phase to LOS displacement
dis = sqrt(gps(:,3).^2+gps(:,5).^2+gps(:,7).^2);
% dis = [gps(:,3);gps(:,5);gps(:,7)];


% Determine subsampling factor for faster plotting
if length(dis) > 400000 && length(dis) < 1000000
    sampling = 2;
elseif length(dis) > 1000000
    sampling = 5;
else
    sampling = 1;
end

% Plot wrapped dataset
figure
quiver(gps(1:sampling:end,1),gps(1:sampling:end,2),gps(1:sampling:end,5),gps(1:sampling:end,3));
axis xy
axis equal
title('GPS')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
% gps=repmat(gps(:,1:2),3,1);
%% Give choice if semi-variogram should be calculated in rectangular area of interest or over the entire image with a mask

choice = questdlg('Would you like to select a rectangular area or mask out a region?', 'Option:', 'Rectangle', 'Mask','Rectangle');

switch choice
    case 'Rectangle'
        dis('Select rectangular area using mouse.')
        bounds = getrect;
        maxLon = bounds(1);
        minLon = bounds(1)+bounds(3);
        minLat = bounds(2);
        maxLat = bounds(2)+bounds(4);
        ixSubset = find(gps(:,2)>minLat & gps(:,2)<maxLat & gps(:,1)<minLon & gps(:,1)>maxLon);
    case 'Mask'
        dis('Draw closed polygon using mouse.')
        polyMask=impoly;
        pos=getPosition(polyMask);
        in = inpolygon(gps(:,1),gps(:,2),pos(:,1),pos(:,2));
        ixSubset = find(in == 0);
end

% Extract subset from rectangular area or after masking
subset = dis(ixSubset);
llon = gps(ixSubset,1);
llat = gps(ixSubset,2);

% Display subregion from selection
% figure('Position', [1, 1, 1200, 1000]);
figure
subplot(2,3,1)
scatter(llon(:),llat(:),[],dis(ixSubset),'.')
colormap(cmapSeismo)
axis xy
axis equal
axis tight
title('Selected region, NON-DETRENDED')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

%% Remove linear trend from subregion
sll = [llon'; llat']';
xy = llh2local(sll',refPoint);
xy = xy*1000;

A = [xy' ones([length(xy) 1])];

coeff = lscov(A,subset);
deramped = subset - A*coeff;

%% Display trend and subregion after removal of trend
subplot(2,3,2)
scatter(llon(:),llat(:),[],A(:,:)*coeff,'.')
colormap(cmapSeismo)
axis xy
axis equal
axis tight
title('Selected region, ESTIMATED TREND')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

subplot(2,3,3)
scatter(llon(:),llat(:),[],deramped(:),'.')
colormap(cmapSeismo)
axis xy
axis equal
axis tight
title('Selected region, DETRENDED')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

%% Calculate and display variogram before plane removal
subplot(2,3,4)
variog = variogram(xy',double(subset),'plotit',true,'subsample',3000);
title('Semi-variogram, NON-DETRENDED')

% Calculate and display variogram after detrending
variogDtrnd = variogram(xy',double(deramped),'plotit',false,'subsample',3000,'nrbins',30);

%% Fit exponential function to experimental variogram and display
subplot(2,3,5)
[a,c,n] = variogramfit(variogDtrnd.distance,variogDtrnd.val,20000,1e-04,variogDtrnd.num, 'model', 'exponential', 'nugget', 1);
title('Semi-variogram and fit, DETRENDED')

h =subplot(2,3,6);
set(h,'visible','off')
text(0.1,1.0,'Fitted exponential semi-variogram parameters:','FontSize',14)
text(0.1,0.8,['Sill:  ', num2str(c)],'FontSize',14)
text(0.1,0.6,['Range:  ', num2str(a)],'FontSize',14)
text(0.1,0.4,['Nugget:  ', num2str(n)],'FontSize',14)

% Print variogram exponential fit parameters to screen
disp(['Sill:  ',num2str(c)])
disp(['Range:  ',num2str(a)])
disp(['Nugget:  ',num2str(n)])

