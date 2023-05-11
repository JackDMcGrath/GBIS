function plotGPSDataModel(gps, geo, invpar, invResults, modelinput, saveName, saveflag, vis)

% Function to generate plot with comparison between GPS data and model
%
% Usage: plotGPSDataModel(gps, geo, invpar, invResults, modelinput, saveName, saveflag)
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Markov chain Monte Carlo algorithm incorporating the Metropolis alghoritm
% (e.g., Mosegaard & Tarantola, JGR,(1995).
%
% by Marco Bagnardi and Andrew Hooper (COMET, University of Leeds)
% Email: M.Bagnardi@leeds.ac.uk
% Reference: Bagnardi and Hooper, 2018
%
% The function uses third party software.
% =========================================================================
% Last update: 11/05/2023
%%

if nargin < 8
    vis='on';
end

global outputDir  % Set global variables

obsGPS = llh2local([gps.ll';zeros(1,length(gps.ll(:,1)))], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
obsGPS = [obsGPS; zeros(1,size(obsGPS,2))]; % Add zeros to third column of observation matrix
nObsGPS = size(obsGPS,2); % Determine number of entries in GPS observation matrix

% Display GPS vectors
scalell = [max(obsGPS(1,:))+5000; min(obsGPS(2,:))-5000; 0]; % add coordinates of legend
hscalebar = max(round(max(abs(gps.displacements(1:2,:)')),3)); % Determine length of scalebar
vscalebar = round(max(abs(gps.displacements(3,:))),3); % Determine length of scalebar

%% Generate plot

if isfield(geo,'faulttracefile')
    fault=dlmread(geo.faulttracefile);
    fault=llh2local([fault, zeros(size(fault(:,1)))]', geo.referencePoint)*1000;
else
    fault=[nan;nan];
end

figure('Position', [1, 1, 1200, 1000],'Visible',vis);
subplot(1,2,1)
hold on
plot(obsGPS(1,:), obsGPS(2,:), 'k.', 'MarkerSize', 10)
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)], [gps.displacements(1,:),-hscalebar], [gps.displacements(2,:),0], 1, ...
'Color','b','LineWidth',1,'MaxHeadSize',0.03);
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)-1e4], [gps.displacements(1,:),-hscalebar], [gps.displacements(2,:),0], 1, ...
    'Color','k','LineWidth',1,'MaxHeadSize',0.03);%,'Marker','s')
axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('GPS horizontal displacements (data:black - model:red)')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
text(scalell(1), scalell(2)-1500,[num2str(hscalebar*1000),' mm'], 'HorizontalAlignment','Right') % Add scalebar
drawnow
hold on
plot(fault(1,:),fault(2,:),'r','LineWidth',2)

% obsGPS(:,end) = []; % remove coordinates of legend
% gps.displacements(:,end) = []; % remove displacements for legend
gps.ve=gps.displacements(:,1);
gps.vn=gps.displacements(:,2);

%% Calculate forward model using optimal source parameters
faults= [find(strcmp(invpar.model,'FAUL')) find(strcmp(invpar.model,'DIKE')) find(strcmp(invpar.model,'SILL')) find(strcmp(invpar.model,'FHIN')) find(strcmp(invpar.model,'SPLT')) find(strcmp(invpar.model,'BSPT')) find(strcmp(invpar.model,'BACK')) find(strcmp(invpar.model,'DHIN'))];
points= [find(strcmp(invpar.model,'MOGI')) find(strcmp(invpar.model,'MCTG')) find(strcmp(invpar.model,'YANG')) find(strcmp(invpar.model,'PENN'))];

modGPS = forwardGPSModel(obsGPS',invpar,invResults,modelinput); % Calculate modelled displacements
modGPSss.ve=modGPS(1,:);
modGPSss.vn=modGPS(2,:);

% % Add back legend
% obsGPS(:,end+1) = [max(obsGPS(1,:))+5000; min(obsGPS(2,:))-5000; 0]; % add coordinates of legend
% gps.displacements(:,end+1) = [-scalebar 0 0]; % add "displacements" of legend

if find(strcmp(invpar.model,'ARCT')) ~= 0
    H_off=invResults.model.mIx(find(strcmp(invpar.model,'ARCT')))+3;
    [gps_dist] = dist2fault_trace(invpar.fault',obsGPS(1:2,:)); % Distances in m
    dist=gps_dist-invResults.model.optimal(H_off);
    scatter(obsGPS(1,:),obsGPS(2,:),20,dist,'filled')
    load('vik.mat')
    [vik]=crop_cmap(vik,[min(dist(:)) max(dist(:))],0);
    colormap(vik)
    colorbar
elseif find(strcmp(invpar.model,'DIST')) ~= 0
    center=invResults.model.optimal(invResults.model.mIx(find(strcmp(invpar.model,'DIST')))+2);
    width=invResults.model.optimal(invResults.model.mIx(find(strcmp(invpar.model,'DIST')))+3);
    dist=invpar.zdist-center;
    inshear=(abs(dist)<=(0.5*width));
    scatter(obsGPS(1,:),obsGPS(2,:),20,dist,'filled')
    scatter(obsGPS(1,inshear),obsGPS(2,inshear),20,dist(inshear),'filled','MarkerEdgeColor','k') % Plot black circle around points within shear zone
    load('vik.mat')
    [vik]=crop_cmap(vik,[min(dist(:)) max(dist(:))],0);
    colormap(vik)
    colorbar
elseif find(strcmp(invpar.model,'DHIN')) ~= 0
    l_edge=invResults.model.mIx(find(strcmp(invpar.model,'DHIN')))+14;
    dist=invpar.zdist-mean([invResults.model.optimal(l_edge),invpar.right]);
    scatter(obsGPS(1,:),obsGPS(2,:),20,dist,'filled')
    load('vik.mat')
    [vik]=crop_cmap(vik,[min(dist(:)) max(dist(:))],0);
    colormap(vik)
    colorbar
end

% Plot modelled displacements (Horizontal)
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)-1e4],[modGPS(1,:),-hscalebar],[modGPS(2,:),0],1,'Color','r','LineWidth',1,'MaxHeadSize',0.03)%,'Marker','s')

if  size(faults,2)>0
    for ii=1:size(faults,2)
        if strcmpi(invpar.model{ii},'FAUL') || strcmpi(invpar.model{ii},'DIKE') || strcmpi(invpar.model{ii},'SILL') || strcmpi(invpar.model{ii},'BACK')
            center=invResults.optimalmodel{faults(ii)}(6:7);
            fault_top=[sind(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1),cosd(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1)];
            fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2))];
            plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center(2),'g','LineWidth',2)
            plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'go','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'FHIN') || strcmpi(invpar.model{ii},'DHIN') || strcmpi(invpar.model{ii},'XHIN')
            PlotHingedFaultsGeometry(invResults.optimalmodel{faults(ii)});
        end
        if strcmpi(invpar.model{ii},'SPLT') % using plotting from splitFault.m
            plotsplit(invResults.optimalmodel{ii});
            plot(invResults.optimalmodel{ii}(6),invResults.optimalmodel{ii}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'BSPT') % using plotting from splitFault.m
            plotsplithinge(invResults.optimalmodel{ii});
        end
    end
end

if  size(points,2)>0
    for ii=1:size(points,2)
        plot(invResults.optimalmodel{points(ii)}(1),invResults.optimalmodel{points(ii)}(2),'b*','MarkerSize',10)
    end
end

subplot(1,2,2)
hold on
plot(obsGPS(1,:), obsGPS(2,:), 'k.', 'MarkerSize', 10)
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)], [gps.displacements(1,:)-modGPS(1,:), -hscalebar], [gps.displacements(2,:)-modGPS(2,:), 0], 1, ...
    'Color','b','LineWidth',1,'MaxHeadSize',0.03);%,'Marker','s')
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)-1e4], [gps.displacements(1,:)-modGPS(1,:), -hscalebar], [gps.displacements(2,:)-modGPS(2,:),0], 1, ...
    'Color','k','LineWidth',1,'MaxHeadSize',0.03);%,'Marker','s')

axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('Horizontal Residuals')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
text(scalell(1), scalell(2)-1500,[num2str(hscalebar*1000),' mm'], 'HorizontalAlignment','Right') % Add scalebar
drawnow
hold on
plot(fault(1,:),fault(2,:),'r','LineWidth',2)

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Data_Model_horizontal'],'-dpng')
end
drawnow

% Plot modelled displacements (Vertical)
figure('Position', [1, 1, 1200, 1000],'Visible',vis);
subplot(1,2,1)
hold on
plot(obsGPS(1,:), obsGPS(2,:), 'k.', 'MarkerSize', 10)
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)], [zeros(size(gps.displacements(1,:))),-vscalebar], [gps.displacements(3,:),0], 1, ...
'Color','b','LineWidth',1,'MaxHeadSize',0.03);
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)-1e4], [zeros(size(gps.displacements(1,:))),-vscalebar], [gps.displacements(3,:),0], 1, ...
    'Color','k','LineWidth',1,'MaxHeadSize',0.03);%,'Marker','s')

axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('GPS vertical displacements (data:black - model:red)')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
text(scalell(1), scalell(2)-1500,[num2str(vscalebar*1000),' mm'], 'HorizontalAlignment','Right') % Add scalebar
drawnow
hold on
plot(fault(1,:),fault(2,:),'r','LineWidth',2)

quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)-1e4], [zeros(size(modGPS(1,:))), -vscalebar],[modGPS(3,:),0],1,'Color','r','LineWidth',1,'MaxHeadSize',0.03)%,'Marker','s')
if  size(faults,2)>0
    for ii=1:size(faults,2)
        if strcmpi(invpar.model{ii},'FAUL') || strcmpi(invpar.model{ii},'DIKE') || strcmpi(invpar.model{ii},'SILL') || strcmpi(invpar.model{ii},'BACK')
            center=invResults.optimalmodel{faults(ii)}(6:7);
            fault_top=[sind(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1),cosd(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1)];
            fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2))];
            plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center(2),'g','LineWidth',2)
            plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'go','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'FHIN')
            PlotHingedFaultsGeometry(invResults.optimalmodel{faults(ii)});
        end
        
        if strcmpi(invpar.model{ii},'SPLT') % using plotting from splitFault.m
            plotsplit(invResults.optimalmodel{ii});
            plot(invResults.optimalmodel{ii}(6),invResults.optimalmodel{ii}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'BSPT') % using plotting from splitFault.m
            plotsplithinge(invResults.optimalmodel{ii});
        end
    end
end

if  size(points,2)>0
    for ii=1:size(points,2)
        plot(invResults.optimalmodel{points(ii)}(1),invResults.optimalmodel{points(ii)}(2),'b*','MarkerSize',10)
    end
end

subplot(1,2,2)
% figure('Visible',vis);
hold on
plot(obsGPS(1,:), obsGPS(2,:), 'k.', 'MarkerSize', 10)
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)], [zeros(size(gps.displacements(1,:))), -vscalebar], [gps.displacements(3,:)-modGPS(3,:),0], 1, ...
    'Color','b','LineWidth',1,'MaxHeadSize',0.03);%,'Marker','s')
quiver([obsGPS(1,:),scalell(1)], [obsGPS(2,:),scalell(2)-1e4], [zeros(size(gps.displacements(1,:))), -vscalebar], [gps.displacements(3,:)-modGPS(3,:),0], 1, ...
    'Color','k','LineWidth',1,'MaxHeadSize',0.03);%,'Marker','s')
axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('Vertical Residuals')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
text(scalell(1), scalell(2)-1500,[num2str(vscalebar*1000),' mm'], 'HorizontalAlignment','Right') % Add scalebar
drawnow
hold on
plot(fault(1,:),fault(2,:),'r','LineWidth',2)


if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Data_Model_vertical'],'-dpng')
end
drawnow

%% Plot contour
plot_contour = 0;
if plot_contour == 1
figure('Position', [1, 1, 1200, 1000],'Visible',vis);
plot(obsGPS(1,:), obsGPS(2,:),'k.');
axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('Modelled GPS vertical displacements (mm/yr)')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);% xlim([-1e4 2e4]),ylim([4.5e4 7.5e4])
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
hold on

if  size(faults,2)>0
    for ii=1:size(faults,2)
        if strcmpi(invpar.model{ii},'FAUL') || strcmpi(invpar.model{ii},'DIKE') || strcmpi(invpar.model{ii},'SILL') || strcmpi(invpar.model{ii},'FHIN') || strcmpi(invpar.model{ii},'BACK')
            center=invResults.optimalmodel{faults(ii)}(6:7);
            fault_top=[sind(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1),cosd(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1)];
            fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2))];
            plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center(2),'g','LineWidth',2)
            plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'go','MarkerSize',10,'MarkerFaceColor','g')
            if strcmpi(invpar.model{ii},'FHIN')
                center2=mean([[fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1)]+center(1);[fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2)]+center(2)]');
                fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10))];
                plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center2(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center2(2),'b--','LineWidth',2)
                plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
            end
        end
        if strcmpi(invpar.model{ii},'SPLT') % using plotting from splitFault.m
            plotsplit(invResults.optimalmodel{ii});
            plot(invResults.optimalmodel{ii}(6),invResults.optimalmodel{ii}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'BSPT') % using plotting from splitFault.m
            plotsplithinge(invResults.optimalmodel{ii});
        end
    end
end

plot(fault(1,:),fault(2,:),'r','LineWidth',2)

%%
rh=max(modGPS(3,:)*1000);

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

lines=floor(min(modGPS(3,:))*1000):b:ceil(max(modGPS(3,:)*1000));


Xmat=[min(obsGPS(1,:)):100:max(obsGPS(1,:))];
Ymat=[min(obsGPS(2,:)):100:max(obsGPS(2,:))];
[Xmat,Ymat]=meshgrid(Xmat,Ymat);
Z=scatteredInterpolant(obsGPS(1,:)',obsGPS(2,:)',(modGPS(3,:)*1000)');
Zmat=Z(Xmat,Ymat);
c=convhull(obsGPS(1,:),obsGPS(2,:));
[in,on]=inpolygon(Xmat,Ymat,obsGPS(1,c),obsGPS(2,c));
out=((in+on)==0);
Zmat(out)=NaN;
[m,C]=contour(Xmat,Ymat,Zmat,lines);
clabel(m,C,lines(1:3:end))
colorbar
drawnow
%%

if  size(points,2)>0
    for ii=1:size(points,2)
        plot(invResults.optimalmodel{points(ii)}(1),invResults.optimalmodel{points(ii)}(2),'b*','MarkerSize',10)
    end
end

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Contour_Model_vertical'],'-dpng')
end
drawnow

%% Plot residual contour

figure('Position', [1, 1, 1200, 1000],'Visible',vis);

axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('Residual East Displacements (Obs-Mod) (mm/yr)')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
hold on


%%
vresid=(gps.displacements(1,:)-modGPS(1,:))*1000;
rh=range(vresid);

if rh>1000
    b=100;
elseif rh>100
    b=10;
elseif rh>50
    b=5;
elseif rh>10
    b=2;
elseif rh>5;
    b=1;
elseif rh>1;
    b=0.5;
else
    b=0.1;
end

lines=round(min(vresid)):b:round(max(vresid));

Xrange=[min(obsGPS(1,:)):100:max(obsGPS(1,:))];
Yrange=[min(obsGPS(2,:)):100:max(obsGPS(2,:))];
[Xmat,Ymat]=meshgrid(Xrange,Yrange);
Z=scatteredInterpolant(obsGPS(1,:)',obsGPS(2,:)',vresid');
Zmat=Z(Xmat,Ymat);
c=convhull(obsGPS(1,:),obsGPS(2,:));
[in,on]=inpolygon(Xmat,Ymat,obsGPS(1,c),obsGPS(2,c));
out=((in+on)==0);
Zmat(out)=NaN;
A=imagesc(Xrange,Yrange,Zmat);
set(A,'AlphaData',~isnan(Zmat));
caxis([min(Zmat(:)) max(Zmat(:))])
colorbar

if  size(faults,2)>0
    for ii=1:size(faults,2)
        if strcmpi(invpar.model{ii},'FAUL') || strcmpi(invpar.model{ii},'DIKE') || strcmpi(invpar.model{ii},'SILL') || strcmpi(invpar.model{ii},'FHIN') || strcmpi(invpar.model{ii},'BACK')
            center=invResults.optimalmodel{faults(ii)}(6:7);
            fault_top=[sind(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1),cosd(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1)];
            fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2))];
            plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center(2),'g','LineWidth',2)
            plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'go','MarkerSize',10,'MarkerFaceColor','g')
            if strcmpi(invpar.model{ii},'FHIN')
                center2=mean([[fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1)]+center(1);[fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2)]+center(2)]');
                fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10))];
                plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center2(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center2(2),'b--','LineWidth',2)
                plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
            end
        end
        if strcmpi(invpar.model{ii},'SPLT') % using plotting from splitFault.m
            plotsplit(invResults.optimalmodel{ii});
            plot(invResults.optimalmodel{ii}(6),invResults.optimalmodel{ii}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'BSPT') % using plotting from splitFault.m
            plotsplithinge(invResults.optimalmodel{ii});
        end
    end
end

plot(fault(1,:),fault(2,:),'r','LineWidth',2)
[m,C]=contour(Xmat,Ymat,Zmat,lines,'-black');
clabel(m,C,lines(1:2:end))
plot(obsGPS(1,:), obsGPS(2,:),'k.');

% load('~/scripts/cpts/vik/vik.mat')
load('vik.mat')
[vik]=crop_cmap(vik,[min(Zmat(:)) max(Zmat(:))],0);
colormap(vik)

if  size(points,2)>0
    for ii=1:size(points,2)
        plot(invResults.optimalmodel{points(ii)}(1),invResults.optimalmodel{points(ii)}(2),'b*','MarkerSize',10)
    end
end

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Contour_Resid_east'],'-dpng')
end

drawnow
%% Plot residual contour

figure('Position', [1, 1, 1200, 1000],'Visible',vis);

axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('Residual North Displacements (Obs-Mod) (mm/yr)')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
hold on


%%
vresid=(gps.displacements(2,:)-modGPS(2,:))*1000;
rh=range(vresid);

if rh>1000
    b=100;
elseif rh>100
    b=10;
elseif rh>50
    b=5;
elseif rh>10
    b=2;
elseif rh>5;
    b=1;
elseif rh>1;
    b=0.5;
else
    b=0.1;
end

lines=round(min(vresid)):b:round(max(vresid));

Xrange=[min(obsGPS(1,:)):100:max(obsGPS(1,:))];
Yrange=[min(obsGPS(2,:)):100:max(obsGPS(2,:))];
[Xmat,Ymat]=meshgrid(Xrange,Yrange);
Z=scatteredInterpolant(obsGPS(1,:)',obsGPS(2,:)',vresid');
Zmat=Z(Xmat,Ymat);
c=convhull(obsGPS(1,:),obsGPS(2,:));
[in,on]=inpolygon(Xmat,Ymat,obsGPS(1,c),obsGPS(2,c));
out=((in+on)==0);
Zmat(out)=NaN;
A=imagesc(Xrange,Yrange,Zmat);
set(A,'AlphaData',~isnan(Zmat));
caxis([min(Zmat(:)) max(Zmat(:))])
colorbar

if  size(faults,2)>0
    for ii=1:size(faults,2)
        if strcmpi(invpar.model{ii},'FAUL') || strcmpi(invpar.model{ii},'DIKE') || strcmpi(invpar.model{ii},'SILL') || strcmpi(invpar.model{ii},'FHIN') || strcmpi(invpar.model{ii},'BACK')
            center=invResults.optimalmodel{faults(ii)}(6:7);
            fault_top=[sind(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1),cosd(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1)];
            fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2))];
            plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center(2),'g','LineWidth',2)
            plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'go','MarkerSize',10,'MarkerFaceColor','g')
            if strcmpi(invpar.model{ii},'FHIN')
                center2=mean([[fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1)]+center(1);[fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2)]+center(2)]');
                fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10))];
                plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center2(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center2(2),'b--','LineWidth',2)
                plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
            end
        end
        if strcmpi(invpar.model{ii},'SPLT') % using plotting from splitFault.m
            plotsplit(invResults.optimalmodel{ii});
            plot(invResults.optimalmodel{ii}(6),invResults.optimalmodel{ii}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'BSPT') % using plotting from splitFault.m
            plotsplithinge(invResults.optimalmodel{ii});
        end
    end
end

plot(fault(1,:),fault(2,:),'r','LineWidth',2)
[m,C]=contour(Xmat,Ymat,Zmat,lines,'-black');
clabel(m,C,lines(1:2:end))
plot(obsGPS(1,:), obsGPS(2,:),'k.');

% load('~/scripts/cpts/vik/vik.mat')
load('vik.mat')
[vik]=crop_cmap(vik,[min(Zmat(:)) max(Zmat(:))],0);
colormap(vik)

if  size(points,2)>0
    for ii=1:size(points,2)
        plot(invResults.optimalmodel{points(ii)}(1),invResults.optimalmodel{points(ii)}(2),'b*','MarkerSize',10)
    end
end

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Contour_Resid_north'],'-dpng')
end

drawnow
%% Plot residual contour

figure('Position', [1, 1, 1200, 1000],'Visible',vis);

axis equal;
ax = gca;
grid on
ax.Layer = 'top';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
xlabel('X distance from local origin (m)')
ylabel('Y distance from local origin (m)')
title('Residual Vertical Displacements (Obs-Mod) (mm/yr)')
xlim([min(obsGPS(1,:))-10000 max(obsGPS(1,:))+10000]);
ylim([min(obsGPS(2,:))-10000 max(obsGPS(2,:))+10000]);
hold on


%%
vresid=(gps.displacements(3,:)-modGPS(3,:))*1000;
rh=range(vresid);

if rh>1000
    b=100;
elseif rh>100
    b=10;
elseif rh>50
    b=5;
elseif rh>10
    b=2;
elseif rh>5;
    b=1;
elseif rh>1;
    b=0.5;
else
    b=0.1;
end

lines=round(min(vresid)):b:round(max(vresid));

Xrange=[min(obsGPS(1,:)):100:max(obsGPS(1,:))];
Yrange=[min(obsGPS(2,:)):100:max(obsGPS(2,:))];
[Xmat,Ymat]=meshgrid(Xrange,Yrange);
Z=scatteredInterpolant(obsGPS(1,:)',obsGPS(2,:)',vresid');
Zmat=Z(Xmat,Ymat);
c=convhull(obsGPS(1,:),obsGPS(2,:));
[in,on]=inpolygon(Xmat,Ymat,obsGPS(1,c),obsGPS(2,c));
out=((in+on)==0);
Zmat(out)=NaN;
A=imagesc(Xrange,Yrange,Zmat);
set(A,'AlphaData',~isnan(Zmat));
caxis([min(Zmat(:)) max(Zmat(:))])
colorbar

if  size(faults,2)>0
    for ii=1:size(faults,2)
        if strcmpi(invpar.model{ii},'FAUL') || strcmpi(invpar.model{ii},'DIKE') || strcmpi(invpar.model{ii},'SILL') || strcmpi(invpar.model{ii},'FHIN') || strcmpi(invpar.model{ii},'BACK')
            center=invResults.optimalmodel{faults(ii)}(6:7);
            fault_top=[sind(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1),cosd(invResults.optimalmodel{faults(ii)}(5))*0.5*invResults.optimalmodel{faults(ii)}(1)];
            fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(4))*invResults.optimalmodel{faults(ii)}(2))];
            plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center(2),'g','LineWidth',2)
            plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'go','MarkerSize',10,'MarkerFaceColor','g')
            if strcmpi(invpar.model{ii},'FHIN')
                center2=mean([[fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1)]+center(1);[fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2)]+center(2)]');
                fault_base=[sind(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10)),cosd(invResults.optimalmodel{faults(ii)}(5)+90)*(cosd(invResults.optimalmodel{faults(ii)}(11))*invResults.optimalmodel{faults(ii)}(10))];
                plot([-fault_top(1), fault_top(1), fault_top(1)-fault_base(1),-fault_top(1)-fault_base(1),-fault_top(1)]+center2(1),[-fault_top(2), fault_top(2),fault_top(2)-fault_base(2),-fault_top(2)-fault_base(2),-fault_top(2)]+center2(2),'b--','LineWidth',2)
                plot(invResults.optimalmodel{faults(ii)}(6),invResults.optimalmodel{faults(ii)}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
            end
        end
        if strcmpi(invpar.model{ii},'SPLT') % using plotting from splitFault.m
            plotsplit(invResults.optimalmodel{ii});
            plot(invResults.optimalmodel{ii}(6),invResults.optimalmodel{ii}(7),'g*','MarkerSize',10,'MarkerFaceColor','g')
        end
        if strcmpi(invpar.model{ii},'BSPT') % using plotting from splitFault.m
            plotsplithinge(invResults.optimalmodel{ii});
        end
    end
end

plot(fault(1,:),fault(2,:),'r','LineWidth',2)
[m,C]=contour(Xmat,Ymat,Zmat,lines,'-black');
clabel(m,C,lines(1:2:end))
plot(obsGPS(1,:), obsGPS(2,:),'k.');

% load('~/scripts/cpts/vik/vik.mat')
load('vik.mat')
[vik]=crop_cmap(vik,[min(Zmat(:)) max(Zmat(:))],0);
colormap(vik)

if  size(points,2)>0
    for ii=1:size(points,2)
        plot(invResults.optimalmodel{points(ii)}(1),invResults.optimalmodel{points(ii)}(2),'b*','MarkerSize',10)
    end
end

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Contour_Resid_vertical'],'-dpng')
end

drawnow


end
%% Plot Pointcloud

figure('Visible',vis);
hold on
title('Observed vs Modelled GPS')
xlabel('GPS Observation (mm)')
ylabel('Modelled GPS (mm)')
plot(gps.displacements(1,:)*1e3,modGPS(1,:)*1e3,'r.')
plot(gps.displacements(2,:)*1e3,modGPS(2,:)*1e3,'b.')
plot(gps.displacements(3,:)*1e3,modGPS(3,:)*1e3,'g.')
legend('East','North','Up','Location','NorthWest','Autoupdate','off')
modrange=[floor(min(modGPS(:)*1e3)) ceil(max(modGPS(:)*1e3))];
obsrange=[floor(min(gps.displacements(:)*1e3)) ceil(max(gps.displacements(:)*1e3))];
try
xlim(obsrange)
catch
    xlim(modrange)
end
try
ylim(modrange)
catch
    ylim([-5 5])
end
plot([-1e3 1e3],[-1e3 1e3],'k--')

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Data_Model_scatter'],'-dpng')
end
drawnow

figure('Visible',vis);
hold on
title('HIstogram of distribution of Verticals')
xlabel('Uplift Rate (mm)')
edges=-10.125:0.25:15.125;
hist(gps.displacements(3,:)*1e3,edges)

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/Vertical_Histogram'],'-dpng')
end
drawnow

%% Plot histogram

vresid=(gps.displacements-modGPS)*1000;
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

E_edges=round(min(vresid(1,:)):b:max(vresid(1,:))); if size(E_edges,2)==1; E_edges(2)=E_edges+b; end
N_edges=round(min(vresid(2,:)):b:max(vresid(2,:))); if size(N_edges,2)==1; N_edges(2)=N_edges+b; end
U_edges=round(min(vresid(3,:)):b:max(vresid(3,:))); if size(U_edges,2)==1; U_edges(2)=U_edges+b; end

histogram(vresid(1,:),E_edges,'FaceColor','r')
histogram(vresid(2,:),N_edges,'FaceColor','b')
histogram(vresid(3,:),U_edges,'FaceColor','g')

legend(sprintf('East- Mean: %.1f, Std: %.1f',mean(vresid(1,:)),std(vresid(1,:))),...
    sprintf('North- Mean: %.1f, Std: %.1f',mean(vresid(2,:)),std(vresid(2,:))),...
    sprintf('Up- Mean: %.1f, Std: %.1f',mean(vresid(3,:)),std(vresid(3,:))),'Location','NorthWest')

xlabel('Residual (mm)')
title('Residual: Observed GPS - Modelled GPS')

if saveflag=='y'
    print(gcf,[outputDir,'/Figures/GPS_Data_Model_hist'],'-dpng')
end

drawnow

