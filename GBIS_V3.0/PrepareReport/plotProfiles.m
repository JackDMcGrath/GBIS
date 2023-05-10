function plotProfiles(geo,GPS,insar,inputFile,invpar,invResults,model,modelInput,obs,nObs,vis,fidHTML,onlyVertFlag)
global outputDir

XHINix=find(ismember(invpar.model,'XHIN'))
Utot=zeros(size(obs));
for i = 1:length(XHINix)
mFunc=invResults.optimalmodel{XHINix(i)};     % Select source model parameters from all
U = hingedFaultsWithOffset(mFunc,obs(1:2,:),modelInput.nu);   % Calculate 3D displacements
Utot=Utot+U;
end

U=Utot;

if sum(ismember(invpar.model,'ARCT')) > 0
mFunc=invResults.optimalmodel{find(ismember(invpar.model,'ARCT'))};
V = arctan(invpar.dist,mFunc(1:4));
U(1,:) = U(1,:) + sind(geo.fault_strike)*V;
U(2,:) = U(2,:) + cosd(geo.fault_strike)*V;
end

U = U*1e3;

paramod=sind(geo.fault_strike)*U(1,:) + cosd(geo.fault_strike)*U(2,:);
perpmod=cosd(geo.fault_strike)*U(1,:) + sind(geo.fault_strike)*U(2,:);

if onlyVertFlag==0
parains=sind(geo.fault_strike)*-insar{1}.dLos + cosd(geo.fault_strike)*-insar{2}.dLos;
perpins=cosd(geo.fault_strike)*-insar{1}.dLos + sind(geo.fault_strike)*-insar{2}.dLos;
end
paragps=sind(geo.fault_strike)*GPS.displacements(1,:) + cosd(geo.fault_strike)*GPS.displacements(2,:);
perpgps=cosd(geo.fault_strike)*GPS.displacements(1,:) + sind(geo.fault_strike)*GPS.displacements(2,:);

if exist('insar') && size(GPS.ll,1)>0
    figure('Position', [1, 1, 3600, 1000],'Visible',vis);
    subplot(1,3,1)
    if onlyVertFlag==0
    plot(invpar.dist(insar{1}.ix),-insar{1}.dLos*1e3,'.')
    end
    hold on
    plot(invpar.dist(GPS.ix),GPS.displacements(1,:)*1e3,'k.','MarkerSize',15)
    plot(invpar.dist,U(1,:),'.')
    title('East')
    
    subplot(1,3,2)
    if onlyVertFlag==0
    plot(invpar.dist(insar{2}.ix),-insar{2}.dLos*1e3,'.')
    end
    hold on
    plot(invpar.dist(GPS.ix),GPS.displacements(2,:)*1e3,'k.','MarkerSize',15)
    plot(invpar.dist,U(2,:),'.')
    title('North')
    
    subplot(1,3,3)
    if onlyVertFlag==0
    plot(invpar.dist(insar{3}.ix),insar{3}.dLos*1e3,'.')
    else
    plot(invpar.dist(insar{1}.ix),insar{1}.dLos*1e3,'.')
    end
    hold on
    plot(invpar.dist(GPS.ix),GPS.displacements(3,:)*1e3,'k.','MarkerSize',15)
    plot(invpar.dist,U(3,:),'.')
    title('Vert')
    
    sgtitle('ENU Profiles')
    
    drawnow

% Save image as png
img = getframe(gcf);
imwrite(img.cdata, [outputDir,'/Figures/ENUProfileModel.png']);

% Add image to html report
fprintf(fidHTML, '%s\r\n', '<hr>');
fprintf(fidHTML, '%s\r\n', '<H3>Convergence plots</H3>');
fprintf(fidHTML, '%s\r\n', ['<img src="Figures/ENUProfileModel.png','" alt="HTML5 Icon">']);
    
    
end

if exist('insar') && size(GPS.ll,1)>0
    figure('Position', [1, 1, 3600, 1000],'Visible',vis);
    subplot(1,3,1)
    if onlyVertFlag==0
    plot(invpar.dist(insar{1}.ix),parains(insar{1}.ix)*1e3,'.')
    end
    hold on
    plot(invpar.dist(GPS.ix),paragps*1e3,'k.','MarkerSize',15)
    plot(invpar.dist,paramod,'.')
    title('Fault Parallel')
    
    subplot(1,3,2)
    if onlyVertFlag==0
    plot(invpar.dist(insar{2}.ix),perpins(insar{1}.ix)*1e3,'.')
    end
    hold on
    plot(invpar.dist(GPS.ix),perpgps*1e3,'k.','MarkerSize',15)
    plot(invpar.dist,perpmod,'.')
    title('Fault Perp')
    
    subplot(1,3,3)
    if onlyVertFlag==0
    plot(invpar.dist(insar{3}.ix),insar{3}.dLos*1e3,'.')
        else
    plot(invpar.dist(insar{1}.ix),insar{1}.dLos*1e3,'.')
    end
    hold on
    plot(invpar.dist(GPS.ix),GPS.displacements(3,:)*1e3,'k.','MarkerSize',15)
    plot(invpar.dist,U(3,:),'.')
    title('Vert')
    
    sgtitle('Fault Parallel Profiles')
    
    drawnow

% Save image as png
img = getframe(gcf);
imwrite(img.cdata, [outputDir,'/Figures/FPProfileModel.png']);

% Add image to html report
fprintf(fidHTML, '%s\r\n', '<hr>');
fprintf(fidHTML, '%s\r\n', '<H3>Convergence plots</H3>');
fprintf(fidHTML, '%s\r\n', ['<img src="Figures/FPProfileModel.png','" alt="HTML5 Icon">']);
    
    
end


