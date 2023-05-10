function [] = fitVariogramGPS(inputFile,step,velcode,subsample)

% Function to fit an exponential function to the data isotropic (semi-)variogram
% Adapted to fit to inversion data by Jack
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
% close all
% clc
if nargin<2
    step=1;
end

if nargin<3
    velcode=['ENU'];
end

if nargin<4
    subsample=3000;
end

% Load dataset
disp('Ingesting data to estimate (semi-)variogram ...')
data = load(inputFile,'inv');
data=data.inv;

% Create colormap
cmapSeismo = colormap_cpt('GMT_seis.cpt', 100);    % GMT 'Seismo' colormap for wrapped data

% Find a local reference point
refPoint = [min(data.lonlat(:,1)), min(data.lonlat(:,2))]; % Determine local reference point
% opts.Interpreter='tex';
% opts.Default='Up';
% choice = questdlg('Which velocity would you like to check?', 'Option:', 'East', 'North', 'Up',opts);
%
% switch choice
%     case 'East'
%         los = data.vx*1e-3;                              % Convert to displacement in m
%     case 'North'
%         los = data.vy*1e-3;                              % Convert to displacement in m
%     case 'All'
%         los = sqrt(((data.vx*1e-3).^2) + ((data.vy*1e-3).^2) + ((data.vz*1e-3).^2));
%     case 'Up'
%         los = data.vz*1e-3;                              % Convert to displacement in m
% end

% Determine subsampling factor for faster plotting
if length(data.lonlat) > 400000 && length(data.lonlat) < 1000000
    sampling = 2;
elseif length(data.lonlat) > 1000000
    sampling = 5;
else
    sampling = 1;
end

% Plot dataset
figure
scatter(data.lonlat(1:sampling:end,1), data.lonlat(1:sampling:end,2), [], data.vz,'.');
colormap(cmapSeismo)
axis xy
axis equal
title('Wrapped Interferogram')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
colorbar

data.lonlat=data.lonlat(1:step:end,:);

%% Give choice if semi-variogram should be calculated in rectangular area of interest or over the entire image with a mask
auto_thresh=1e11;
if subsample >= auto_thresh
    fprintf('At least %.0f subsamples selected. Automatically doing 50-100km.\nSod off and leave me to it. Come back later.\n', auto_thresh)
    choice='50-100km';
else
    choice = questdlg('Would you like to select a rectangular area or mask out a region?', 'Option:', 'Rectangle', 'Mask','50-100km', 'Rectangle');
end

switch choice
    case 'Rectangle'
        choice = questdlg('Rectangle by area or closed polygon?', 'Option:', 'Area', 'Polygon', 'Rectangle');
        switch choice
            case 'Area'
                disp('Select rectangular area using mouse.')
                bounds = getrect;
                maxLon = bounds(1);
                minLon = bounds(1)+bounds(3);
                minLat = bounds(2);
                maxLat = bounds(2)+bounds(4);
                ixSubset = find(data.lonlat(:,2)>minLat & data.lonlat(:,2)<maxLat & data.lonlat(:,1)<minLon & data.lonlat(:,1)>maxLon);
            case 'Polygon'
                choice = questdlg('Default polygon?', 'Option:', 'Default', 'Manual', 'Rectangle');
                    switch choice
                        case 'Default'
                            pos=[399,5106;446,5148;494,5089;445,5049];
                        case 'Manual'
                            disp('Draw closed polygon using mouse.')
                            polyMask=impoly;
                            pos=getPosition(polyMask);
                    end
                in = inpolygon(data.lonlat(:,1),data.lonlat(:,2),pos(:,1),pos(:,2));
                ixSubset = find(in==1);
        end
    case 'Mask'
        disp('Draw closed polygon using mouse.')
        polyMask=impoly;
        pos=getPosition(polyMask);
        in = inpolygon(data.lonlat(:,1),data.lonlat(:,2),pos(:,1),pos(:,2));
        ixSubset = find(in == 0);
    case '50-100km'
        disp('Using predefined region')
        ixSubset = find(inpolygon(data.lonlat(:,1),data.lonlat(:,2),[376,433,523,510],[5097,5076,5139,5189]));
        
end

% Extract subset from rectangular area or after masking
llon = data.lonlat(ixSubset,1);
llat = data.lonlat(ixSubset,2);
c=convhull(llon,llat);
hold on
plot(llon(c),llat(c));
drawnow

for ii=1:length(velcode)
    if strcmpi(velcode(ii),'E')
        los = data.vx(1:step:end)*1e-3;
        direction='East';
    elseif strcmpi(velcode(ii),'N')
        los = data.vy(1:step:end)*1e-3;
        direction='North';
    elseif strcmpi(velcode(ii),'U')
        los = data.vz(1:step:end)*1e-3;
        direction='Up';
    end
    
    % Extract subset from rectangular area or after masking
    subset = los(ixSubset);
    
    % Display subregion from selection
    % figure('Position', [1, 1, 1200, 1000]);
    figure('Name',direction);
    subplot(2,3,1)
    scatter(llon(:),llat(:),[],subset,'.')
    colormap(cmapSeismo)
    axis xy
    axis equal
    axis tight
    title('Selected region, NON-DETRENDED')
    xlabel('Longitude (UTM km)')
    ylabel('Latitude (UTM km)')
    colorbar
    
    %% Remove linear trend from subregion
    xy = data.lonlat(ixSubset,:)';
    xy = xy*1000;
    
    A = [xy' ones([length(xy) 1])];
    
    coeff = lscov(A,subset');
    
    if sum(isnan(coeff)) ~=0
        coeff([1,2,3])=0;
        fprintf('No Trend Found!\n')
    end
    
    deramped = subset' - A*coeff;
    %% Display trend and subregion after removal of trend
    subplot(2,3,2)
    scatter(llon(:),llat(:),[],A(:,:)*coeff,'.')
    colormap(cmapSeismo)
    axis xy
    axis equal
    axis tight
    title('Selected region, ESTIMATED TREND')
    xlabel('Longitude (UTM km)')
    ylabel('Latitude (UTM km)')
    colorbar
    
    subplot(2,3,3)
    scatter(llon(:),llat(:),[],deramped(:),'.')
    colormap(cmapSeismo)
    axis xy
    axis equal
    axis tight
    title('Selected region')
    xlabel('Longitude (UTM km)')
    ylabel('Latitude (UTM km)')
    colorbar
    
    %% Calculate and display variogram before plane removal
    subplot(2,3,4)
    variog = variogram(xy',double(subset)','plotit',true,'subsample',subsample);
    title('Semi-variogram, NON-DETRENDED')
    
    % Calculate and display variogram after detrending
    variogDtrnd = variogram(xy',double(deramped),'plotit',false,'subsample',subsample,'nrbins',30);
    
    %% Fit exponential function to experimental variogram and display
    subplot(2,3,5)
    [a,c,n] = variogramfit(variogDtrnd.distance,variogDtrnd.val,20000,1e-04,variogDtrnd.num, 'model', 'exponential', 'nugget', 1);
    title('Semi-variogram and fit')
    
    h =subplot(2,3,6);
    set(h,'visible','off')
    text(0.1,1.0,'Fitted exponential semi-variogram parameters:','FontSize',14)
    text(0.1,0.8,['Sill:  ', num2str(c)],'FontSize',14)
    text(0.1,0.6,['Range:  ', num2str(a)],'FontSize',14)
    text(0.1,0.4,['Nugget:  ', num2str(n)],'FontSize',14)
    
    % Print variogram exponential fit parameters to screen
    disp(['Results for ',direction,' velocity with stepsize ',num2str(step)])
    disp(['Sill:  ',num2str(c)])
    disp(['Range:  ',num2str(a)])
    disp(['Nugget:  ',num2str(n)])
end
