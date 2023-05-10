function [insar, obs, nObs] = loadInsarData(insar, geo, cmap, plotfig)

% Function to ingest and subsample InSAR data from pre-prepared *.mat file
%
% Usage: [insar, obs, nObs] = loadInsarData(insar, geo, cmap)
% Input Parameters:
%       insar: structure with insar data and related information
%       geo: structure with local coordinates origin and bounding box
%       cmap: colormaps for plotting
%
% Output Parameters:
%       insar: structure with added subsampled data vector and radar look
%       parameters
%       obs: coordinates of observation points after subsampling
%       nObs: number of observation points
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

if nargin<4
    plotfig=1;
end

%% Initialise variables
global outputDir  % Set global variables

nPoints = 0;
LonLat = zeros(0,2);

%% Start loading data

if plotfig==0
    vis='off';
    fprintf('InSAR figure not plotted. Change flag to 1\n')
elseif plotfig==1
    vis='on';
end

for i = 1:length(insar)
    loadedData = load(insar{i}.dataPath); % load *.mat file
    
    % Apply bounding box and remove data points outside the AOI
    iOutBox = find(loadedData.Lon < geo.boundingBox(1) | loadedData.Lon > geo.boundingBox(3) | loadedData.Lat > geo.boundingBox(2) | loadedData.Lat < geo.boundingBox(4));
    if sum(iOutBox)>0
        loadedData.Phase(iOutBox) = [];
        loadedData.Lat(iOutBox) = [];
        loadedData.Lon(iOutBox) = [];
        loadedData.Heading(iOutBox) = [];
        loadedData.Incidence(iOutBox) = [];
    end
    
    convertedPhase = (loadedData.Phase / (4*pi) )  * insar{i}.wavelength;    % Convert phase from radians to m
    los = single(-convertedPhase);  % Convert to Line-of-sight displacement in m
    ll = [single(loadedData.Lon) single(loadedData.Lat)];   % Create Longitude and Latitude 2Xn matrix
    xy = llh2local(ll', geo.referencePoint);    % Transform from geografic to local coordinates
    
    nPointsThis = size(ll, 1);   % Calculate length of current InSAR data vector
    xy = double([(1:nPointsThis)', xy'*1000]);   % Add ID number column to xy matrix with local coordinates
    
    % Extract filename to be included in figure names
    [path, name, ext] = fileparts(insar{i}.dataPath);
    
    % Plot wrapped interferogram
    figure('Position', [1, 1, 700, 700],'Visible',vis);
    plotInsarWrapped(xy, los, insar{i}.wavelength, cmap, name);
    saveas(gcf, [outputDir, '/Figures/Wrapped_', name, '.png'])
    
    % Plot unwrapped interferogram
    figure('Position', [1, 1, 700, 700],'Visible',vis);
    plotInsarUnwrapped(xy, los, cmap, name);
    saveas(gcf,[outputDir,'/Figures/Unwrapped_',name,'.png'])
    
    %% Run data vector subsampling using Quadtree and display
    if insar{i}.quadtreeThresh~=0
        [nb, err, nPts, centers, dLos, polys, xLims, yLims] = quadtree(xy, los', insar{i}.quadtreeThresh, 1000, 1,vis); % Run Quadtree on los vector
        c = max(abs([min(los), max(los)])); % Calculate maximu value for symmetric colormap
        caxis([-c c])
        axis equal; axis tight;
        cbar = colorbar; ylabel(cbar, 'Line-of-sight displacement m','FontSize', 14);
        colormap(cmap.redToBlue)
        xlabel('X distance from local origin (m)','FontSize', 14)
        ylabel('Y distance from local origin (m)','FontSize', 14)
        t = title(['Subsampled data. Number of data points used:', num2str(nb)],'FontSize', 18);
        set(t,'Position',get(t,'Position')+[0 1000 0]);
        drawnow
        saveas(gcf, [outputDir,'/Figures/Subsampled_', name, '.png'])
        
        %% Extract radar look vector information for subsampled points
        disp 'Extracting radar look vector parameters ...'
        [dHeading] = quadtreeExtract(xy, loadedData.Heading, xLims(:,1:4), yLims(:,1:4));    % Extract heading angle values based on quadtree partition
        disp(['Mean heading angle: ', num2str(mean(dHeading)), ' degrees'])
        [dIncidence] = quadtreeExtract(xy, loadedData.Incidence, xLims(:,1:4), yLims(:,1:4));  % Extract values based on quadtree partition
        disp(['Mean incidence angle: ', num2str(mean(dIncidence)), ' degrees'])
        disp(['Max and min LoS displacement in m:', num2str(max(dLos)), '  ', num2str(min(dLos))])
        disp 'Extracting observation points height ...'
        
        pts.xy = [(1:nb)', centers];    % Generate Nx3 [# x y] matrix with local coordinates of data points (centers of Quadtree cells)
        pts.LonLat = local2llh(pts.xy(:,2:3)'/1000, geo.referencePoint)'; % Convert x and y coordinates into Lon Lat coordinates
        
        % Include Quadtree results into insar structure
        insar{i}.obs = pts.xy(:,2:3);
        insar{i}.dLos = dLos;
        insar{i}.dHeading = dHeading;
        insar{i}.dIncidence = dIncidence;
    else
        fprintf('Not Applying Quadtree subsampling\n')
        insar{i}.obs = xy(:,2:3);
        pts.LonLat = local2llh(xy(:,2:3)'/1000, geo.referencePoint)'; % Convert x and y coordinates into Lon Lat coordinates
        insar{i}.dLos = los';
        insar{i}.dHeading = loadedData.Heading';
        insar{i}.dIncidence = loadedData.Incidence';
        nb=length(insar{i}.dLos);
    end
    
    %% Create inverse of covariance matrix
%     obs = pts.xy(:,2:3);
    obs = insar{i}.obs;    

    [X1,X2] = meshgrid(obs(:,1)); % Create square matrices of Xs
    [Y1,Y2] = meshgrid(obs(:,2)); % Create square matrices of Ys
    H = sqrt((X1-X2).^2 + (Y1 - Y2).^2); % Calculate distance between points
    
    % Assign default values if sill, range and nugget are not provided
    if  ~isfield(insar{i},'sillExp');
        disp 'sillExp value not found, assigning default value 0.04^2'
        insar{i}.sillExp = 0.04^2;
    end
    
    if ~isfield(insar{i},'nugget');
        disp 'nugget value not found, assigning default value 0.002^2'
        insar{i}.nugget = 0.002^2;
    end
    
    if ~isfield(insar{i},'range');
        disp 'range value not found, assigning default value 0.04^2'
        insar{i}.range = 5000;
    end
    
    covarianceMatrix = insar{i}.sillExp * exp(-3*H/insar{i}.range) + insar{i}.nugget*eye(nb); % Calculate covariance matrix for exponential model with nugget
    insar{i}.invCov = inv(covarianceMatrix); % Calculate inverse of covariance matrix
    clear X1 X2 Y1 Y2 obs H covarianceMatrix
    %%
    insar{i}.ix = nPoints+1:nPoints+nb; % Extract index of data points in obs vector for this interferogram
    nPoints = nPoints + nb; % Number of data points
    LonLat = [LonLat; pts.LonLat]; % Longitude Latitude matrix
end

obs = llh2local([LonLat'; zeros(1,nPoints)],geo.referencePoint')*1000;
obs = [obs; zeros(1, size(obs,2))]; % Coordinates of observation points
nObs = size(obs,2); % Total number of observation points

