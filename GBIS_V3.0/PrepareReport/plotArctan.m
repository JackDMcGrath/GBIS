function plotArctan(geo,GPS,insar,inputFile,invpar,invResults,model,modelInput,obs,nObs,vis,fidHTML,onlyVertFlag)
global outputDir

if exist('insar') && size(GPS.ll,1)==0
    figure('Position', [1, 1, 1200, 1000],'Visible',vis);
    subplot(1,2,1)
    if onlyVertFlag==0
    para=sind(geo.fault_strike)*-insar{1}.dLos + cosd(geo.fault_strike)*-insar{2}.dLos;
    quiver(obs(1,1:nObs/size(insar,2)),obs(2,1:nObs/size(insar,2)),sind(geo.fault_strike)*para,cosd(geo.fault_strike)*para)
    end
    fault=dlmread(geo.faulttracefile);
    fault=llh2local([fault, zeros(size(fault(:,1)))]', geo.referencePoint)*1000;
    hold on
    plot(fault(1,:),fault(2,:))
    xlim([min(obs(1,:))-5e3 max(obs(1,:))+5e3])
    ylim([min(obs(2,:))-5e3 max(obs(2,:))+5e3])
    axis equal
    
    subplot(1,2,2)
    plot(invpar.dist(1:nObs/size(insar,2)),para,'.')
    hold on
    mFunc=invResults.optimalmodel{find(ismember(invpar.model,'ARCT'))}; % Select source model parameters from all; Opening set to 0;
    V = arctan(invpar.dist,mFunc(1:4));
    plot(invpar.dist,V,'.')
    sgtitle('Fault Parallel Component, and Shear Zone Contribution')
    
elseif exist('insar') && size(GPS.ll,1)>0
    if onlyVertFlag==0
    para=sind(geo.fault_strike)*-insar{1}.dLos + cosd(geo.fault_strike)*-insar{2}.dLos;
    end
    paragps=sind(geo.fault_strike)*GPS.displacements(1,:) + cosd(geo.fault_strike)*GPS.displacements(2,:);
    figure('Position', [1, 1, 1200, 1000],'Visible',vis);
    subplot(1,2,1)
    if onlyVertFlag==0
    quiver(obs(1,insar{1}.ix),obs(2,insar{1}.ix),sind(geo.fault_strike)*para,cosd(geo.fault_strike)*para)
    end
    hold on
    plot(obs(1,end-size(GPS.ll,1)+1:end),obs(2,end-size(GPS.ll,1)+1:end),'k.','MarkerSize',10)
    fault=dlmread(geo.faulttracefile);
    fault=llh2local([fault, zeros(size(fault(:,1)))]', geo.referencePoint)*1000;
    plot(fault(1,:),fault(2,:))
    xlim([min(obs(1,:))-5e3 max(obs(1,:))+5e3])
    ylim([min(obs(2,:))-5e3 max(obs(2,:))+5e3])
    axis equal
    
    subplot(1,2,2)
    if onlyVertFlag==0
    plot(invpar.dist(insar{1}.ix),para,'.')
    end
    hold on
    plot(invpar.dist(end-size(GPS.ll,1)+1:end),paragps,'k.','MarkerSize',15)
    mFunc=invResults.optimalmodel{find(ismember(invpar.model,'ARCT'))}; % Select source model parameters from all; Opening set to 0;
    V = arctan(invpar.dist,mFunc(1:4));
    plot(invpar.dist,V,'.')
    sgtitle('Fault Parallel Component, and Shear Zone Contribution')
else
    paragps=sind(geo.fault_strike)*GPS.displacements(1,:) + cosd(geo.fault_strike)*GPS.displacements(2,:);
    figure('Position', [1, 1, 1200, 1000],'Visible',vis);
    plot(invpar.dist,paragps,'k.','MarkerSize',15)
    hold on
    mFunc=invResults.optimalmodel{find(ismember(invpar.model,'ARCT'))}; % Select source model parameters from all; Opening set to 0;
    V = arctan(invpar.dist,mFunc(1:4));
    plot(invpar.dist,V,'.')
    title('Fault Parallel Component, and Shear Zone Contribution')
    
end

drawnow

% Save image as png
img = getframe(gcf);
imwrite(img.cdata, [outputDir,'/Figures/arctanModel.png']);

% Add image to html report
fprintf(fidHTML, '%s\r\n', '<hr>');
fprintf(fidHTML, '%s\r\n', '<H3>Convergence plots</H3>');
fprintf(fidHTML, '%s\r\n', ['<img src="Figures/arctanModel.png','" alt="HTML5 Icon">']);
