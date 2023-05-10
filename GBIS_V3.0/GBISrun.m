function [] = GBISrun(inputFileName, insarDataCode, gpsDataCode, modelCode, nRuns, skipSimulatedAnnealing,displayfigs,finalReport)

%%  Geodetic Bayesian Inversion Software (GBIS)
%   Software for the Bayesian inversion of geodetic data.
%   Markov chain Monte Carlo algorithm incorporating the Metropolis algorithm
%
%%  =========================================================================
%   Usage: GBISrun(inputFileName, insarDataCode, gpsDataFlag, modelCode, nRuns, skipSimulatedAnnealing)
%
%   inputFileName:  name and extension of input file (e.g., 'VolcanoExercise.inp')
%
%   insarDataCode:  select data to use (e.g., [1,3] to use InSAR data
%                   with insarID = 1 and 3; insarID specified in input file *.inp).
%                   Leave empty (e.g.,[]) if no InSAR data.
%
%   gpsDataFlag:    'y' to use GPS data, 'n' to not use GPS data
%
%   modelCode:      select forward models to use;
%                   'M' for Mogi source (point source, [Mogi, 1958])
%                   'T' for McTigue source (finite sphere, [McTigue, 1987])
%                   'Y' for Yang source (prolate spheroid, [Yang et al., 1988])
%                   'S' for rectangular horizontal sill with uniform opening [Okada, 1985]
%                   'P' for penny-shaped crack [Fialko et al., 2001]
%                   'D' for rectangular dipping dike with uniform opening [Okada, 1985]
%                   'F' for rectangular dipping fault with uniform slip [Okada, 1985]
%
%                   Custom made sources:
%                   'H' for two hinged rectangular dikes (hinged at the bottom edge of the upper dike)
%                   'B' for two hinged rectangular faults (hinged at the bottom edge of the upper fault)
%                   'W' for a torn fault, with two sections of different dip (tear location on the surface trace not locking depth)
%                   'C' for a split hinged fault
%                   'N' for a backlip fault
%
%   nRuns:          number of iterations (samples) to be performed (e.g., 1000000)
%
%   skipSimulatedAnnealing: 'y' to skip initial Simulated Annealing phase
%
%   Example: >> GBISrun('VolcanoExercise.inp',[1,2],'y','MD',1e06, 'n')
%   =========================================================================
%   Contacts
%   by Marco Bagnardi and Andrew Hooper (COMET, University of Leeds)
%   Email: M.Bagnardi@leeds.ac.uk
%   Reference: TBA (Bagnardi and Hooper, in prep.)
%   Last update: 12/05/2017

%% Check number of input arguments and return error if not sufficient

if nargin == 0
    help GBISrun;
    return;
end

if nargin < 6
    disp('#############################################')
    disp('####### Not enough input arguments. #########')
    disp('#############################################')
    disp(' ')
    disp('Type "help GBISrun" for more information');
    disp(' ')
    return;
end

if nargin < 7
    displayfigs=1;
end

if nargin < 8
    finalReport=1;
end

if length(gpsDataCode) > 0;
    gpsDataFlag = 'y';
else
    gpsDataFlag = 'n';
end

%% Start timer and initialise global variables
%clc % Clean screen
tic; % Start timer

clear  outputDir  % Clear global variables
global outputDir % Set global variables

%% Diplay header
disp('**********************************************')
disp('Geodetic Bayesian Inversion Software (GBIS)')
disp('Software for the Bayesian inversion of geodetic data.')
disp('Markov chain Monte Carlo algorithm incorporating the Metropolis algorithm')
disp(' ')
disp('by Marco Bagnardi and Andrew Hooper (COMET, University of Leeds)')
disp('Email: M.Bagnardi@leeds.ac.uk')
disp('Last update: 12/05/2017')
disp('**********************************************')
disp(' ')
disp(' ')

%% Read input data from input text file *.inp

inputFileID = fopen(inputFileName, 'r');
textLine = fgetl(inputFileID); 

while ischar(textLine)
    eval(textLine)
    textLine = fgetl(inputFileID);
end

fclose(inputFileID);

%% Create output directories and output file name

[inputFile.path, inputFile.name, inputFile.ext] = fileparts(inputFileName); % Extract name
outputDir = inputFile.name; % Name output directory as input file name
    
% if ~isfield(gps,'enu_weight') % Flag to set equal GPS ENU weights if not set in input
%         gps.enu_weight=[1;1;1];
% end

   % InSAR + GPS
if gpsDataFlag == 'y' && ~isempty(insarDataCode)
    
    % Add InSAR dataset ID
    for i = 1:length(insarDataCode)
        insarDataNames{i} = num2str(insarDataCode(i));
    end
    
    % Add GPS
      for i = 1:length(gpsDataCode)
        gpsDataNames{i} = num2str(gpsDataCode(i));
    end
    saveName = ['invert_',strjoin(insarDataNames, '_'),'_GPS',strjoin(gpsDataNames, '_')];
    
    % Add GPS weight
%     saveName = sprintf('%s%.0f%.0f%.0f',saveName,gps.enu_weight(1),gps.enu_weight(2),gps.enu_weight(3));
    
    % Add models
    for i = 1:length(modelCode)
        saveName = [saveName,'_',modelCode(i)];
    end
    
    if exist('out_extension') && size(out_extension,2)>0
        saveName = [saveName,'_',out_extension];
    end
    
    % Create output directories
    outputDir = [outputDir,'/',saveName];
    disp(['Output directory: ', outputDir])
    mkdir(outputDir)
    mkdir([outputDir,'/Figures']) % Create directory for Figures
    copyfile(inputFileName,outputDir) % Copy input file to output directory
    
    % Add .mat extension
    saveName = [saveName,'.mat'];

    % InSAR only
elseif gpsDataFlag == 'n' && ~isempty(insarDataCode)
    
    % Add InSAR dataset ID
    for i = 1:length(insarDataCode)
        insarDataNames{i} = num2str(insarDataCode(i));
    end

    saveName = ['invert_',strjoin(insarDataNames, '_')];
    
    % Add models
    for i = 1:length(modelCode)
        saveName = [saveName,'_',modelCode(i)];
    end
    
    if exist('out_extension') && size(out_extension,2)>0
        saveName = [saveName,'_',out_extension];
    end
    
    % Create output directories
    outputDir = [outputDir,'/',saveName];
    disp(['Output directory: ', outputDir])
    mkdir(outputDir)
    mkdir([outputDir,'/Figures']) % Create directory for Figures
    copyfile(inputFileName,outputDir) % Copy input file to output directory
    
    % Add .mat extension
    saveName = [saveName,'.mat'];
    
    % GPS only
elseif gpsDataFlag == 'y' && isempty(insarDataCode)
    
    % Add GPS
          for i = 1:length(gpsDataCode)
        gpsDataNames{i} = num2str(gpsDataCode(i));
    end
    saveName = ['invert_GPS',strjoin(gpsDataNames, '_')];        
    
    % Add models
    for i = 1:length(modelCode)
        saveName = [saveName,'_',modelCode(i)];
    end
    
    if exist('out_extension') && size(out_extension,2)>0
        saveName = [saveName,'_',out_extension];
    end
            
    % Create output directories
    outputDir = [outputDir,'/',saveName];
    disp(['Output directory: ', outputDir])
    mkdir(outputDir)
    mkdir([outputDir,'/Figures']) % Create directory for Figures
    copyfile(inputFileName,outputDir) % Copy input file to output directory
    
    
    % Add .mat extension
    saveName = [saveName,'.mat'];
end

%% Initialise variables

nObs = 0;   % Initialise number of observations variable
obs  = [];  % Initialise observation points (x,y,z) matrix

%% Ingest InSAR data

% Create colormaps for plotting InSAR data (call third party colormap_cpt.m function and GMT *.cpt files)
cmap.seismo = colormap_cpt('/nfs/see-fs-02_users/eejdm/scripts/insar/GBIS/GBIS_V1.0/GBIS/ThirdPartyScripts/GMT_seis.cpt', 100);    % GMT 'Seismo' colormap for wrapped data
cmap.redToBlue = colormap_cpt('/nfs/see-fs-02_users/eejdm/scripts/insar/GBIS/GBIS_V1.0/GBIS/ThirdPartyScripts/polar.cpt', 100);    % Red to Blue colormap for unwrapped data

% Select InSAR datasets to use from list in input file
if ~isempty(insarDataCode)
    disp('InSAR datasets used in this inversion:')
    
    for i = 1:length(insarDataCode)
        selectedInsarData{i} = insar{insarDataCode(i)};
        disp(insar{insarDataCode(i)}.dataPath)    % display filename of InSAR dataset
    end
    
    % Create subsampled InSAR datasets for inversion
    disp(' ')
    disp('Ingesting InSAR data and performing Quadtree subsampling ...')
    [insar, obsInsar, nObsInsar] = loadInsarData(selectedInsarData, geo, cmap,displayfigs); % Load and subsample InSAR data
    nObs = nObsInsar;   % Add number of InSAR data points to total number of points
    obs  = obsInsar;    % Add InSAR observation points to observation points
else
    disp(' ')
    disp 'No InSAR datasets will be used in this inversion.'
    insar = [];
end

%% Ingest GPS data
GPS = [];
nGPS = 0;
if gpsDataFlag == 'y'
    for i = 1:length(gpsDataCode)
    disp(' ')
    fprintf('GPS dataset %.0f used in this inversion:\n',gpsDataCode(i))
    disp(gps{gpsDataCode(i)}.dataPath) % display filename of GPS data file
    disp(' ')
    disp('Ingesting GPS data ...')
    [GPS{i}, obsGps, nObsGps] = loadGpsData(gps{gpsDataCode(i)}, geo,displayfigs);
    disp([num2str(nObsGps), ' GPS sites will be used in this inversion'])
    GPS{i}.ix = [nObs+1:nObs+nObsGps]; % Index of GPS data in observation vector
    nObs = [nObs + nObsGps];   % Add number of GPS data points to total number of points
    obs = [obs, obsGps];       % Add GPS observation points to observation points
    nGPS = nGPS+nObsGps;
    end
    gps = GPS;
    
    if length(gpsDataCode) > 1
        if displayfigs==1
            vis='on';
        else
            vis='off';
        end
        gpsobs=[];
        gpsdisp=[];
        for i = 1:length(gpsDataCode)
            gpsobs=[gpsobs,obs(:,gps{i}.ix)];
            gpsdisp=[gpsdisp,gps{i}.displacements];
        end
        scalebar = abs(round(max(gpsdisp(:))/3,3));
        
        
        figure('Visible',vis)
        subplot(1,2,1)
        quiver(gpsobs(1,:), gpsobs(2,:), gpsdisp(1,:), gpsdisp(2,:), 1, 'LineWidth', 1, 'MaxHeadSize', 0.03, 'Marker', 's')
        axis equal;
        ax = gca;
        grid on
        ax.Layer = 'top';
        ax.Box = 'on';
        ax.LineWidth = 1.5;
        ax.GridLineStyle = '--';
        xlabel('X distance from local origin (m)'); ylabel('Y distance from local origin (m)')
        title('GPS horizontal displacements')
        xlim([min(gpsobs(1,:))-10000 max(gpsobs(1,:))+10000]);
        ylim([min(gpsobs(2,:))-10000 max(gpsobs(2,:))+10000]);
        text(gpsobs(1,end),gpsobs(2,end)-2000,[num2str(scalebar*1000),' mm'])
        
        subplot(1,2,2)
        quiver(gpsobs(1,:), gpsobs(2,:), zeros(1,size(gpsdisp,2)), gpsdisp(3,:), 1, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 0.03)
        axis equal;
        ax = gca;
        grid on
        ax.Layer = 'top';
        ax.Box = 'on';
        ax.LineWidth = 1.5;
        ax.GridLineStyle = '--';
        xlabel('X distance from local origin (m)'); ylabel('Y distance from local origin (m)')
        title('GPS vertical displacements')
        xlim([min(gpsobs(1,:))-10000 max(gpsobs(1,:))+10000]);
        ylim([min(gpsobs(2,:))-10000 max(gpsobs(2,:))+10000]);
        text(gpsobs(1,end),gpsobs(2,end)-2000,[num2str(scalebar*1000),' mm'])
        
        if strcmpi(vis,'on');drawnow;end
        if ~isempty(outputDir)
            saveas(gcf,[outputDir,'/Figures/GPS_displacements_hor_vert.png'])
        end
        
    end
else
    disp(' ')
    disp 'No GPS datasets will be used in this inversion.'
    gps = [];
end

%% Plug-in here any further type of data to ingest (i.e., differential DEMs)

%% Setup inversion parameters

% Pause and press key to continue
disp ' '
disp '#################   Press any key to continue ...'
% pause

% Define inversion parameters
disp 'Preparing for inversion ...'

invpar.nSave = 1000;    % Save output to file every 1000 iterations (every 10,000 after 20,000 iterations)
invpar.sensitivitySchedule = [1:100:10000,11000:1000:30000,40000:10000:nRuns]; % sensitivity schedule (when to change step sizes)

if skipSimulatedAnnealing == 'y'
    invpar.TSchedule = 1; % No temperature schedule if Simulated Annealing is not performed
else
    invpar.TSchedule = 10.^(3:-0.2:0); % Cooling schedule for Simulated Annealing
end
invpar.TRuns = 1000; % Number of runs at each temperature (Simulated Annealing only)

invpar.nModels = length(modelCode); % number of models used (e.g., 'MD' = 2x models)

invpar.nRuns = nRuns;

if isfield(geo,'fault_strike')
    invpar.fault_strike = geo.fault_strike;
end

% Switch model code to full model name
% ADD HERE ANY NEW CUSTOMISED MODEL

for i = 1:invpar.nModels
    modelName = modelCode(i);
    switch modelName
        case 'M'
            invpar.model{i}='MOGI';    % Mogi source
        case 'T'
            invpar.model{i}='MCTG';    % McTigue source
        case 'P'
            invpar.model{i}='PENN';   % Penny-shaped crack (Fialko 2001)
        case 'Y'
            invpar.model{i}='YANG';    % Yang source
        case 'S'
            invpar.model{i}='SILL';    % Sill simulated as horizontal dislocation (Okada)
        case 'D'
            invpar.model{i}='DIKE';    % Dipping dike dislocation (Okada)
        case 'F'
            invpar.model{i}='FAUL';    % Dipping fault dislocation (Okada)
            
            % Custom made models
            % New Deformation Types
        case 'A'
            invpar.model{i}='ARCT';   % Arctan function
            invpar.fault=readmatrix(geo.faulttracefile); % Load in fault trace to measure distance fault from
            invpar.fault=llh2local([invpar.fault(:,[1 2]), zeros(size(invpar.fault(:,1)))]', geo.referencePoint)*1000;
            [invpar.dist] = dist2fault_trace(invpar.fault',obs([1 2],:)); % Distances in m
        case 'Z'
            invpar.model{i}='DIST';   % Distributed Shear Zone function
            invpar.fault=readmatrix(geo.faulttracefile); % Load in fault trace to measure distance fault from
            invpar.fault=llh2local([invpar.fault(:,[1 2]), zeros(size(invpar.fault(:,1)))]', geo.referencePoint)*1000;
            [invpar.zdist] = dist2fault_trace(invpar.fault',obs([1 2],:))/1000; % Distances in km
        case 'N'
            invpar.model{i}='BACK';   % Backslip
            
            % Combination Models
        case 'H'
            invpar.model{i}='HING';    % Two dikes hinged along L at depth (2x Okada)
        case 'B'
            invpar.model{i}='FHIN';    % Two faults hinged along L at depth (2x Okada)
        case 'O'
            invpar.model{i}='XHIN';    % Two faults hinged along L at depth (2x Okada) with X-offet
        case 'W'
            invpar.model{i}='SPLT';    % One fault torn into 2 differently dipping sections (2x Okada)
        case 'C'
            invpar.model{i}='BSPT';    % One FHIN torn into 2 differently dipping sections (4x Okada)
        case 'X'
            invpar.model{i}='DHIN';    % One FHIN with DIST attatched to the hingepoint at the right edge (2x Okada, 1x Dist)
            invpar.fault=readmatrix(geo.faulttracefile); % Load in fault trace to measure distance fault from
            invpar.fault=llh2local([invpar.fault(:,[1 2]), zeros(size(invpar.fault(:,1)))]', geo.referencePoint)*1000;
            [invpar.zdist] = dist2fault_trace(invpar.fault',obs([1 2],:))/1000; % Distances in km
            % End of custom made models
            
        otherwise
            error 'Invalid model code.'
    end
end

model = prepareModel(modelInput, invpar, insar, gps);

%% Run inversion

[invResults,invpar] = runInversion(geo, gps, insar, invpar, model, modelInput, obs, nObs);

%% Create *.mat file with final results
% disp('=========================================================')
disp('Saving Results')
disp(['Time since start (HH:MM:SS):  ',datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')])
if gpsDataFlag == 'y' && ~isempty(insar)

    cd(outputDir)
    save(saveName, 'insarDataCode', 'geo', 'inputFile', 'invpar', 'gps', 'insar', 'model', 'modelInput', 'invResults', 'obs', 'nObs', 'saveName','-v7.3')
    delete temporary.mat
    cd ../..
    
elseif gpsDataFlag == 'n' && ~isempty(insar)

    cd(outputDir)
    save(saveName, 'insarDataCode', 'geo', 'inputFile', 'invpar', 'insar', 'model', 'modelInput', 'invResults', 'obs', 'nObs', 'saveName','-v7.3')
    delete temporary.mat
    cd ../..
    
elseif gpsDataFlag == 'y' && isempty(insar)

    cd(outputDir)
    save(saveName, 'geo', 'inputFile', 'invpar', 'gps', 'model', 'modelInput', 'invResults', 'obs', 'nObs', 'saveName','-v7.3')
    delete temporary.mat
    cd ../..
end

%% Display inversion duration
disp('=========================================================')
disp('GBIS inversion Complete')
disp(['Time since start (HH:MM:SS):  ',datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')])

if finalReport==1
cd(outputDir)
if nRuns == 1
    generateFinalReport(saveName,0,0)
elseif nRuns <= 1e4
    generateFinalReport(saveName,1,0)
elseif nRuns < 5e5 && nRuns > 1e4
    generateFinalReport(saveName,1e4,0)
else
    generateFinalReport(saveName,nRuns*0.2,0)
end
disp(['Time since start (HH:MM:SS):  ',datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')])
cd ../..
end
