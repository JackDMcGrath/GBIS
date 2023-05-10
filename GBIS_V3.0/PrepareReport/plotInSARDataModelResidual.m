function plotInSARDataModelResidual(insar, geo, invpar, invResults, modelinput, saveName, fidHTML, saveflag,vis)

% Function to generate plot with comparison between InSAR data, model, and
% residuals
%
% Usage: plotInSARDataModelResidual(insar, geo, invpar, invResults, modelinput, saveName, fidHTML, saveflag)
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
% Last update: 09/05/2017
%%

if nargin < 9
    vis='on';
end

global outputDir  % Set global variables

% Create colormaps for plotting InSAR data
cmap.seismo = colormap_cpt('GMT_seis.cpt', 100); % GMT Seismo colormap for wrapped interferograms
cmap.redToBlue = colormap_cpt('polar.cpt', 100); % Red to Blue colormap for unwrapped interferograms

for i=1:length(insar)
    % Load and display DATA
    loadedData = load(insar{i}.dataPath); % load *.mat file
    
    % Apply bounding box and remove data points outside it
    iOutBox = find(loadedData.Lon<geo.boundingBox(1) | loadedData.Lon>geo.boundingBox(3) | loadedData.Lat>geo.boundingBox(2) | loadedData.Lat<geo.boundingBox(4));
    if sum(iOutBox)>0
        loadedData.Phase(iOutBox) = [];
        loadedData.Lat(iOutBox) = [];
        loadedData.Lon(iOutBox) = [];
        loadedData.Heading(iOutBox) = [];
        loadedData.Incidence(iOutBox) = [];
    end
    
    convertedPhase = (loadedData.Phase/(4*pi))*insar{i}.wavelength;    % Convert phase from radians to m
    los = -convertedPhase;  % Convert phase from cm to Line-of-sigth displacement in m
    if loadedData.Incidence(1)==90 % So west and south are positive in scatter plots
        plot_data{i}.los = -los';
    else
        plot_data{i}.los = los';
    end
    Heading = loadedData.Heading;
    Inc = loadedData.Incidence;
    ll = [single(loadedData.Lon) single(loadedData.Lat)];   % Create Longitude and Latitude matrix
    xy = llh2local(ll', geo.referencePoint);    % Transform from geographic to local coordinates
    
    nPointsThis = size(ll,1);   % Calculate length of current data vector
    xy = double([[1:nPointsThis]',xy'*1000]);   % Add ID number column to xy matrix with local coordinates
    
    % Patch scattered data for faster plotting
    edge = round(min(abs(diff(xy(:,3)))))+2; % Size of patch set to minumum distance between points
    if edge < 50
        edge = 50;
    end
    xs = [xy(:,2)'; xy(:,2)'+edge; xy(:,2)'+edge; xy(:,2)']; % Coordinates of four vertex of patch
    ys = [xy(:,3)'; xy(:,3)'; xy(:,3)'+edge; xy(:,3)'+edge];
    
    % Calculate MODEL
    constOffset = 0;
    xRamp = 0;
    yRamp = 0;
    
    if i == 1
        if insar{i}.constOffset == 'y'
            constOffset = invResults.model.mIx(end);
            invResults.model.mIx(end) = invResults.model.mIx(end)+1;
        end
        if insar{i}.rampFlag == 'y'
            xRamp = invResults.model.mIx(end);
            yRamp = invResults.model.mIx(end)+1;
            invResults.model.mIx(end) = invResults.model.mIx(end)+2;
        end
    end
    
    if i > 1
        if insar{i}.constOffset == 'y'
            constOffset = invResults.model.mIx(end);
            invResults.model.mIx(end) = invResults.model.mIx(end)+1;
        end
        if insar{i}.rampFlag == 'y'
            xRamp = invResults.model.mIx(end);
            yRamp = invResults.model.mIx(end)+1;
            invResults.model.mIx(end) = invResults.model.mIx(end)+2;
        end
    end
    
    modLos = forwardInsarModel(insar{i},xy,invpar,invResults,modelinput,geo,Heading,Inc,constOffset,xRamp,yRamp); % Modeled InSAR displacements
    if loadedData.Incidence(1)==90 % So west and south are positive in scatter plots
        plot_data{i}.modLos = -modLos;
    else
        plot_data{i}.modLos = modLos;
    end
    
    % Extract filename to be included in figure name
    [pathstr,name,ext] = fileparts(insar{i}.dataPath);
    load('vik.mat')
    fprintf('Wavelength / 1000. Change if using actual InSAR\n')
    [vikUw.redToBlue]=crop_cmap(vik,[min(min(los),min(modLos))  max(max(los),max(modLos))],0);
    % Display wrapped DATA interferogram at 5.6 cm wavelength
    figure('Position', [1, 1, 1200, 1000],'Visible',vis);
    ax1 = subplot(2,3,1);
    plotInsarWrapped(xy,los, insar{i}.wavelength*1e-3, cmap, 'DATA');
    freezeColors
    
    % Display DATA unwrapped interferogram
    ax2 = subplot(2,3,4);
    plotInsarUnwrapped(xy,los, vikUw, 'DATA');
    c = max(abs([min(los), max(los)])); % Calculate maximum value for symmetric colormap
    %     caxis([-c c])
    caxis([min(min(los),min(modLos))  max(max(los),max(modLos))])
    freezeColors
    
    % Display MODEL wrapped interferogram at 5.6 cm wavelength
    ax3=subplot(2,3,2);
    plotInsarWrapped(xy,modLos',insar{i}.wavelength*1e-3,  cmap, 'MODEL');
    colormap(ax3,cmap.seismo)
    freezeColors
    
    % Display MODEL unwrapped interferogram
    ax4=subplot(2,3,5);
    plotInsarUnwrapped(xy,modLos', vikUw, 'MODEL');
    %     caxis([-c c])
    caxis([min(min(los),min(modLos))  max(max(los),max(modLos))])
    freezeColors
    
    % Display RESIDUAL wrapped interferogram at 5.6 cm wavelength
    residual = los-modLos';
    ax5=subplot(2,3,3);
    plotInsarWrapped(xy,residual, insar{i}.wavelength*1e-3, cmap, 'RESIDUAL');
    freezeColors
    
    % Display RESIDUAL unwrapped interferogram
    ax6=subplot(2,3,6);
    plotInsarUnwrapped(xy,residual, cmap, 'RESIDUAL');
    caxis([-c c])
    freezeColors
    
    colormap(ax1,cmap.seismo)
    colormap(ax2,vikUw.redToBlue)
    colormap(ax3,cmap.seismo)
    colormap(ax4,vikUw.redToBlue)
    colormap(ax5,cmap.seismo)
    colormap(ax6,cmap.redToBlue)
    img = getframe(gcf);
    if saveflag=='y'
        imwrite(img.cdata,[outputDir,'/Figures/InSAR_Data_Model_Residual_',name,'.png']);
        
        % Add image to html report
        fprintf(fidHTML, '%s\r\n', '<BR></BR><H3>Comparison InSAR Data - Model - Residual</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/InSAR_Data_Model_Residual_',name,'.png','" alt="HTML5 Icon">']);
    end
    
    
    insar_ix = i;
end

%% Plot Pointcloud

figure('Visible',vis);
hold on
title('Observed vs Modelled InSAR')
xlabel('InSAR Observation (mm)')
ylabel('Modelled InSAR (mm)')
modrange=[];
obsrange=[];

for i=1:length(insar)
    n=split(insar{i}.dataPath,'/');
    lab{i}=replace(n{end}(1:end-4),'_',' ');
    plot(plot_data{i}.los*1e3,plot_data{i}.modLos*1e3,'.')
    obsrange=[obsrange,plot_data{i}.los];
    modrange=[modrange,plot_data{i}.modLos];
end
legend(lab,'Location','NorthWest','AutoUpdate','Off')
modrange=[floor(min(modrange*1e3)) ceil(max(modrange*1e3))];
obsrange=[floor(min(obsrange*1e3)) ceil(max(obsrange*1e3))];
xlim(obsrange)
try
    ylim(modrange)
catch
    ylim([-5 5])
end
plot([-1e3 1e3],[-1e3 1e3],'k--')

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/InSAR_Data_Model_scatter'],'-dpng')
    % Add image to html report
    fprintf(fidHTML, '%s\r\n', '<BR></BR><H3>Observed vs Modelled InSAR Data</H3>');
    fprintf(fidHTML, '%s\r\n', ['<img src="Figures/InSAR_Data_Model_scatter.png','" alt="HTML5 Icon">']);
end
drawnow

for i=1:length(insar)
    if insar{i}.dIncidence(1)==0
        figure('Visible',vis);
        hold on
        title('Histogram of distribution of Verticals')
        xlabel('Uplift Rate (mm)')
        edges=-10.125:0.25:15.125;
        hist(plot_data{i}.los*1e3,edges)
        
        if saveflag=='y'
            print(gcf,[outputDir,'/Figures/Vertical_InSAR_Histogram'],'-dpng')
            % Add image to html report
            fprintf(fidHTML, '%s\r\n', '<BR></BR><H3>Distribution of Initial Verticals</H3>');
            fprintf(fidHTML, '%s\r\n', ['<img src="Figures/Vertical_InSAR_Histogram.png','" alt="HTML5 Icon">']);
        end
        drawnow
    end
end
%% Plot histogram
obs=[];
mod=[];
try
    for i=1:length(insar)
        obs=[obs;plot_data{i}.los];
        mod=[mod;plot_data{i}.modLos];
    end
    vresid=(obs-mod)*1000;
    rh=range(vresid(:));
    
    if rh>1000
        b=100;
    elseif rh>100
        b=10;
    elseif rh>50
        b=5;
    elseif rh>10
        b=1;
    elseif rh>5;
        b=0.5;
    elseif rh>1;
        b=0.25;
    else
        b=0.1;
    end
    
    figure('Visible',vis);
    hold on
    
    label=[];
    for i=1:length(insar)
        if insar{i}.dIncidence(1)==90
            if insar{i}.dHeading(1)==0
                E_edges=round(min(vresid(i,:)):b:max(vresid(i,:))); if size(E_edges,2)==1; E_edges(2)=E_edges+b; end
                histogram(vresid(i,:),E_edges,'FaceColor','r')
                label{i}=sprintf('East- Mean: %.1f, Std: %.1f',mean(vresid(i,:)),std(vresid(i,:)));
            elseif insar{i}.dHeading(1)==270
                N_edges=round(min(vresid(i,:)):b:max(vresid(i,:))); if size(N_edges,2)==1; N_edges(2)=N_edges+b; end
                histogram(vresid(i,:),N_edges,'FaceColor','b')
                label{i}=sprintf('North- Mean: %.1f, Std: %.1f',mean(vresid(i,:)),std(vresid(i,:)));
            end
        elseif insar{i}.dIncidence(1)==0
            U_edges=round(min(vresid(i,:)):b:max(vresid(i,:))); if size(U_edges,2)==1; U_edges(2)=U_edges+b; end
            histogram(vresid(i,:),U_edges,'FaceColor','g')
            label{i}=sprintf('Up- Mean: %.1f, Std: %.1f',mean(vresid(i,:)),std(vresid(i,:)));
        end
    end
    legend(label,'Location','NorthWest')
    
    xlabel('Residual (mm)')
    title('Residual: Observed GPS - Modelled GPS')
    
    if saveflag=='y'
        print(gcf,[outputDir,'/Figures/InSAR_Data_Model_hist'],'-dpng')
        % Add image to html report
        fprintf(fidHTML, '%s\r\n', '<BR></BR><H3>InSAR Model Residuals Histogram</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/InSAR_Data_Model_hist.png','" alt="HTML5 Icon">']);
    end
catch
    fprintf('Not Printing Histograms\n')
end

drawnow

