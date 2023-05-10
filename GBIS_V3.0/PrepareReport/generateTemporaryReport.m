function generateTemporaryReport(invResFile,burning,insarDataCode)

if nargin < 2
    disp('!!!!!!! Not enough input parameters. !!!!!!!!!')
    return;
end

clear  outputDir  % Clear global variables
global outputDir  % Set global variables

load(invResFile);

invResults.mKeep = mKeep;
invResults.model = model;
invResults.optimalmodel = model.funcOpt;

% outputDir = 'temporary';
% mkdir(outputDir);
% mkdir([outputDir,'/Figures']);

% Create colormaps for plotting InSAR data
cmap.Seismo = colormap_cpt('GMT_seis.cpt', 100); % GMT Seismo colormap for wrapped interferograms
cmap.RnB = colormap_cpt('polar.cpt', 100); % Red to Blue colormap for unwrapped interferograms

nParam = length(invResults.mKeep(:,1));

% Print results
format shortG
disp('    Parameter n.   Optimal         Mean       Median         2.5%        97.5%')
for i = 1:nParam-1
    disp([i invResults.model.optimal(i) mean(invResults.mKeep(i,burning:end-1e04)) median(invResults.mKeep(i,burning:end-1e04)) ...
        prctile(invResults.mKeep(i,burning:end-1e04),2.5) prctile(invResults.mKeep(i,burning:end-1e04),97.5)])
end

% Plot convergence
figure('Position', [1, 1, 1200, 1000]);
for i = 1:nParam-1
    subplot(round(nParam/3),3,i)
    plot([1:500:length(invResults.mKeep(1,:))-10000],invResults.mKeep(i,1:500:end-10000),'r.')
end
drawnow


% Plot histograms and best fitting model
figure('Position', [1, 1, 1200, 1000]);
for i = 1:nParam-1
        subplot(round(nParam/3),3,i)
    h=histogram(invResults.mKeep(i,burning:end-1e04),100,'EdgeColor','none','Normalization','count');
    hold on
    topLim = max(h.Values)+10000;
    plot([invResults.model.optimal(i),invResults.model.optimal(i)],[0,topLim],'r-')
    ylim([0 topLim])
%     xlim([mean(invResults.mKeep(i,burning:end-1e04))-2*std(invResults.mKeep(i,burning:end-1e04)), ...
%         mean(invResults.mKeep(i,burning:end-1e04))+2*std(invResults.mKeep(i,burning:end-1e04))]);
end
drawnow
% % Plot joint probabilities
% figure
% plotmatrix_lower(invResults.mKeep(1:nParam-1,burning:end-1e04)','contour');

choice = questdlg(['Do you want to compare DATA MODEL and RESIDUAL?'], 'Plot?', 'Yes', 'No','Yes');
switch choice
    case 'Yes'
        % Plot GPS data, model
        gps = [];
        
        if ~isempty(gps)
            plotGPSDataModel(gps,geo,invpar, invResults, modelinput, 'y')
        end
        
        % Plot InSAR data, model, residual
        if ~isempty(insarDataCode)
            plotInSARDataModelResidual(insar, geo, invpar, invResults, modelInput, 'y')
        end
    case 'No'
        return
end








