function generateFinalReport(invResFile,burning,plotfigs)

% Function to generate summary text file, results plots, and html report
%
% Usage: generateFinalReport(invResFile,burning)
% Input parameters:
%           invResFile: path and name to file with final results of
%                       inversion
%                       (e.g.,'VolcanoExercise/invert_1__MOGI_DIKE.mat')
%           burning:    number of iterations to ignore in pdf histogram plot
%                       and in computation of mean/median/confidence interval
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
if nargin < 2
    disp('!!!!!!! Not enough input parameters. !!!!!!!!!')
    return;
end

if nargin < 3
    plotfigs=1;
end

if plotfigs==1; vis='on';elseif plotfigs==0; vis='off';end

clear  outputDir  % Clear global variables
global outputDir  % Set global variables

warning('off','all')

load(invResFile);   % load results


%% Because I just cant be bothered putting all the GPS into forloops for plotting
GPS.ll = [];
GPS.displacements = [];
GPS.ix = [];
if exist('gps')
    for i = 1:length(gps)
        GPS.ll = [GPS.ll; gps{i}.ll];
        GPS.displacements = [GPS.displacements, gps{i}.displacements];
        GPS.ix = [GPS.ix, gps{i}.ix];
    end
    gps = GPS;
end


outputDir = pwd;

% Create colormaps for plotting InSAR data
cmap.Seismo = colormap_cpt('GMT_seis.cpt', 100); % GMT Seismo colormap for wrapped interferograms
cmap.RnB = colormap_cpt('polar.cpt', 100); % Red to Blue colormap for unwrapped interferograms

nParam = length(invResults.mKeep(:,1)); % Number of model parameters

% Number of empty cells at the end of mKeep and pKeep
if invpar.nRuns < 10000
    blankCells = 999;
else
    blankCells = 9999;
end

% Check to ensure that all blank cells are identified - for some reason it missed the last row, made the joint probs a mess
% (This test will be a disaster if everything converges at 0)
while sum(invResults.mKeep(:,end-blankCells))==0
    blankCells=blankCells+1;
end


%% Print results to text file
% Print optimal, mean, median, 2.5th and 97.5th percentiles
format shortG

txtName = ['/summary',saveName(7:end-4),'_',num2str(burning),'.txt']; % name of text file
fileID = fopen([outputDir,txtName],'w');
fprintf(fileID,'%s\r\n','GBIS');
fprintf(fileID,'%s\r\n',['Summary for ',saveName(1:end-4)]);
if invpar.nRuns > 1 % Added switch in case only 1 run required (ie for forward models)
    fprintf(fileID,'%s\r\n',['Number of iterations: ',num2str(invpar.nRuns)]);
    fprintf(fileID,'%s\r\n',['Model Rejection Rate: ',num2str(invResults.reject_percent),'%']);
    fprintf(fileID,'%s\r\n',['Burning time (n. of iterations): ',num2str(burning)]);
    fprintf(fileID,'%s\r\n','=============================================================================================================');
    fprintf(fileID,'%s\r\n','MODEL PARAM.      OPTIMAL         MEAN            Median          2.5%            97.5%            Starting');
    
    for i = 1:nParam-1
        fprintf(fileID,'%12s\t %8g\t %8g\t %8g\t %8g\t %8g\t %8g\r\n', ...
            char(invResults.model.parName(i)), invResults.model.optimal(i), ...
            mean(invResults.mKeep(i, burning:end-blankCells)), median(invResults.mKeep(i, burning:end-blankCells)), ...
            prctile(invResults.mKeep(i, burning:end-blankCells),2.5), prctile(invResults.mKeep(i, burning:end-blankCells),97.5),...
            model.m(i));
    end
else
    fprintf(fileID,'%s\r\n',['Forward Model Parameters']);
    fprintf(fileID,'%s\r\n','==================================');
    fprintf(fileID,'%s\r\n','MODEL PARAM.      ASSIGNED');
    
    for i = 1:nParam-1
        fprintf(fileID,'%12s\t %8g\t\n', ...
            char(invResults.model.parName(i)), invResults.model.optimal(i));
    end
end

% Display text on screen
%clc
type([outputDir,txtName])

%% Create report html file
htmlName = ['/report',saveName(7:end-4),'_',num2str(burning),'.html']; % HTML file name

% Print header and summary to file
fidHTML = fopen([outputDir,htmlName],'w');
fprintf(fidHTML, '%s\r\n', '<!DOCTYPE html>');
fprintf(fidHTML, '%s\r\n', '<html>');
fprintf(fidHTML, '%s\r\n', '<head>');
fprintf(fidHTML, '%s\r\n', ['<H1>GBIS Final report for <i>', saveName(1:end-4),'</i></H1>']);
fprintf(fidHTML, '%s\r\n', ['<H3>Results file: <i>',outputDir,'/',saveName,'</i></H3>']);
fprintf(fidHTML, '%s\r\n', ['<p>Number of iterations: ', num2str(invpar.nRuns),'</p>']);
fprintf(fidHTML, '%s\r\n', ['Model Rejection Rate: ',num2str(invResults.reject_percent),'%']);
fprintf(fidHTML, '%s\r\n', ['<p>Burning time (n. of iterations from start): ', num2str(burning),'</p>']);
fprintf(fidHTML, '%s\r\n', '<hr>');
fprintf(fidHTML, '%s\r\n', '<H3>Model parameters</H3>');
fprintf(fidHTML, '%s\r\n', '<style>');
fprintf(fidHTML, '%s\r\n', 'table {font-family: arial, sans-serif; border-collapse: collapse; width:100%%;}');
fprintf(fidHTML, '%s\r\n', 'td, th {border: 1px solid #dddddd;text-align: right;padding: 8px;}');
fprintf(fidHTML, '%s\r\n', 'tr:nth-child(even) {background-color: #bbb;}');
fprintf(fidHTML, '%s\r\n', '</style>');
fprintf(fidHTML, '%s\r\n', '</head>');
fprintf(fidHTML, '%s\r\n', '<body>');
fprintf(fidHTML, '%s\r\n', '<table>');

if invpar.nRuns > 1 % Added switch in case only 1 run required (ie for forward models)
    fprintf(fidHTML, '%s\r\n', '<tr> <th>Parameter</th> <th>Optimal</th> <th>Mean</th> <th>Median</th> <th>2.5%</th> <th>97.5%</th> <th>Starting</th></tr>');
    for i = 1:nParam-1
        fprintf(fidHTML, '%s\r\n',['<tr> <td>',char(invResults.model.parName(i)),...
            '</td> <td>', num2str(invResults.model.optimal(i), '%.4f'), ...
            '</td> <td>', num2str(mean(invResults.mKeep(i,burning:end-blankCells)), '%.3f'), ...
            '</td> <td>', num2str(median(invResults.mKeep(i,burning:end-blankCells)), '%.3f'),...
            '</td> <td>', num2str(prctile(invResults.mKeep(i,burning:end-blankCells),2.5), '%.3f'), ...
            '</td> <td>', num2str(prctile(invResults.mKeep(i,burning:end-blankCells),97.5), '%.3f'), ...
            '</td> <td>', num2str(model.m(i), '%.3f') ...
            '</td> </tr>']);
    end
else
    fprintf(fidHTML, '%s\r\n', '<tr> <th>Parameter</th> <th>Assigned</th>');
    for i = 1:nParam-1
        fprintf(fidHTML, '%s\r\n',['<tr> <td>',char(invResults.model.parName(i)),...
            '</td> <td>', num2str(invResults.model.optimal(i), '%.2f'), ...
            '</td> </tr>']);
    end
end

fprintf(fidHTML, '%s\r\n', '</table> </body>');

%% Plot GPS motions

 if exist('gps')
    % Add image to html report
    fprintf(fidHTML, '%s\r\n', '<hr>');
    fprintf(fidHTML, '%s\r\n', '<H3>GPS Observations</H3>');
    fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_displacements_hor_vert.png','" alt="HTML5 Icon">']);
 end
 
 if ~isempty(strfind(saveName,'invert_1_2_3'))
     onlyVertFlag=0;
 else
     onlyVertFlag=1;
 end
 
 
 if ~isempty(find(ismember(invpar.model,'ARCT')))
     % Just plotting for arctan testing WILL FAIL IF YOU'RE NOT DOING ENU
     % INSAR
     
     plotArctan(geo,GPS,insar,inputFile,invpar,invResults,model,modelInput,obs,nObs,vis,fidHTML,onlyVertFlag)
 end
     %plotProfiles(geo,GPS,insar,inputFile,invpar,invResults,model,modelInput,obs,nObs,vis,fidHTML,onlyVertFlag)
%% Plot convergence of all parameters

if invpar.nRuns > 1 % Added switch in case only 1 run required (ie for forward models)
    
    variable_param=find(model.step); % Doesn't plot fixed parameters (don't know what the final parameter is though?)
    figure('Position', [1, 1, 1200, 1000],'Visible',vis);% xlim([-1e4 2e4]),ylim([4.5e4 7.5e4])
    if invpar.nRuns <= 1e4
        step = 10;
    else
        step = 100;
    end
    
    for i = 1:numel(variable_param)-1
        subplot(ceil(numel(variable_param)/3),3,i)    % Determine poistion in subplot
        plot(1:step:length(invResults.mKeep(1,:))-blankCells, invResults.mKeep(variable_param(i),1:step:end-blankCells),'r.') % Plot one point every 100 iterations
        hold on
        try
            xline(burning,'--'); % xline requires matlab/2018b
        end
        title(invResults.model.parName(variable_param(i)))
    end
    if strcmp(invpar.model{1},'FHIN') || strcmp(invpar.model{1},'XHIN')
        m=invResults.model.optimal;
        subplot(ceil(numel(variable_param)/3),3,numel(variable_param))    % Determine poistion in subplot
        if length(invpar.model) > 1;
            if  strcmp(invpar.model{2},'DIST')
                plot(polyshape([m(16)-0.5*m(17),m(16)-0.5*m(17),m(16)+0.5*m(17),m(16)+0.5*m(17)],[-m(15),-100,-100,-m(15)]),'FaceColor','r')
            elseif strcmp(invpar.model{2},'ARCT')
                if strcmp(invpar.model{1},'FHIN')
                plot([m(17),m(17)]*1e-3,[-m(15)*1e-3,-100],'r')
                elseif strcmp(invpar.model{1},'XHIN')
                    plot([m(18),m(18)]*1e-3,[-m(16)*1e-3,-100],'r')
                end
            end
        end
        hold on
        if strcmp(invpar.model{1},'FHIN')
        plot([m(3)/tand(-m(4)),(m(3)/tand(-m(4)))+(cosd(m(4))*m(2)),120e3]*1e-3,[-m(3),-(m(3)-sind(m(4))*m(2)),-(m(3)-sind(m(4))*m(2))]*1e-3,'b')
        else
            plot([m(3)/tand(-m(4)),m(3)/tand(-m(4))+m(14)]*1e-3,[-m(3),-m(3)]*1e-3,'g')
            plot([m(3)/tand(-m(4))+m(14),(m(3)/tand(-m(4)))+(cosd(m(4))*m(2))+m(14),120e3]*1e-3,[-m(3),-(m(3)-sind(m(4))*m(2)),-(m(3)-sind(m(4))*m(2))]*1e-3,'b')
        end
        title('Geometry')
        axis equal
        ylim([-50 0])
        xlim([0 100])
        hold off
    end
    
    % Save image as png
    img = getframe(gcf);
    imwrite(img.cdata, [outputDir,'/Figures/Convergence_',num2str(burning),'.png']);
    
    % Add image to html report
    fprintf(fidHTML, '%s\r\n', '<hr>');
    fprintf(fidHTML, '%s\r\n', '<H3>Convergence plots</H3>');
    fprintf(fidHTML, '%s\r\n', ['<img src="Figures/Convergence_',num2str(burning),'.png','" alt="HTML5 Icon">']);
    
    %% Plot histograms and optimal values
    
    figure('Position', [1, 1, 1200, 1000],'Visible',vis);
    
    for i = 1:numel(variable_param)-1
        subplot(round(numel(variable_param)/3),3,i) % Determine poistion in subplot
        hold on
        xMin = mean(invResults.mKeep(variable_param(i),burning:end-blankCells))-4*std(invResults.mKeep(variable_param(i),burning:end-blankCells));
        xMax = mean(invResults.mKeep(variable_param(i),burning:end-blankCells))+4*std(invResults.mKeep(variable_param(i),burning:end-blankCells));
        bins = xMin: (xMax-xMin)/50: xMax;
        try
            h = histogram(invResults.mKeep(variable_param(i),burning:end-blankCells),bins,'EdgeColor','none','Normalization','count');
            topLim = max(h.Values);
            plot([invResults.model.optimal(variable_param(i)),invResults.model.optimal(variable_param(i))],[0,topLim+10000],'r-') % Plot optimal value
            ylim([0 topLim+10000])
        catch
            fprintf('Skipped Histogram for %s\n',invResults.model.parName{i}) % Don't plot histograms for values forced to be constant
        end
        hold on
        title(invResults.model.parName(variable_param(i)))
    end
    % Save image as png
    img = getframe(gcf);
    imwrite(img.cdata,[outputDir,'/Figures/PDFs_',num2str(burning),'.png']);
    
    % Add image to html report
    fprintf(fidHTML, '%s\r\n', '<hr>');
    fprintf(fidHTML, '%s\r\n', '<BR></BR><H3>Model parameters posterior probabilities and optimal values</H3>');
    fprintf(fidHTML, '%s\r\n', ['<img src="Figures/PDFs_',num2str(burning),'.png','" alt="HTML5 Icon">']);
    
end

%% Plot joint probabilities
if plotfigs==1
    if invpar.nRuns > 1 % Added switch in case only 1 run required (ie for forward models)
        % Optional
        choice = questdlg('Do you want to plot the joint probabilities?', 'Plot?', 'Yes', 'No','Yes');
        switch choice
            case 'Yes'
                
                figure('Position', [1, 1, 1200, 1000]);
                [~,~,~,P,pax]=plotmatrix_lower(invResults.mKeep(variable_param(1:end-1),burning:end-blankCells)','contour',model.parName(variable_param(1:end-1)));
                for ii=1:numel(variable_param)-1
                    tmp_ax=get(pax(ii),'Xlabel');
                    tmp_ax.String=model.parName(variable_param(ii));
                end
                img1 = getframe(gcf);
                imwrite(img1.cdata,[outputDir,'/Figures/JointProbabilities_',num2str(burning),'.png']);
                
                
                % Add image to html report
                fprintf(fidHTML, '%s\r\n', '<hr>');
                fprintf(fidHTML, '%s\r\n', '<H3>Joint probabilities</H3>');
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/JointProbabilities_',num2str(burning),'.png','" alt="HTML5 Icon">']);
                
            case 'No'
        end
    end
    %% Plot comparison betweem data, model, and residual
    
    % Optional
    choice = questdlg('Do you want to compare DATA MODEL and RESIDUAL?', 'Plot?', 'Yes', 'No','Yes');
    switch choice
        case 'Yes'
            % Plot GPS data, model
            
            if exist('gps')
                plotGPSDataModel(gps,geo,invpar, invResults, modelInput, saveName, 'y')
                
                % Add image to html report
                fprintf(fidHTML, '%s\r\n', '<hr>');
                fprintf(fidHTML, '%s\r\n', '<H3>Comparison GPS Data vs. Model</H3>');
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Data_Model_horizontal.png','" alt="HTML5 Icon">']);
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Data_Model_vertical.png','" alt="HTML5 Icon">']);
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Resid_east.png','" alt="HTML5 Icon">']);
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Resid_north.png','" alt="HTML5 Icon">']);
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Model_vertical.png','" alt="HTML5 Icon">']);
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Resid_vertical.png','" alt="HTML5 Icon">']);
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Data_Model_scatter.png','" alt="HTML5 Icon">']);
                fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Data_Model_hist.png','" alt="HTML5 Icon">']);
            end
            
            % Plot InSAR data, model, residual
            if exist('insarDataCode')
                plotInSARDataModelResidual(insar, geo, invpar, invResults, modelInput, saveName, fidHTML, 'y')
            end
        case 'No'
            return
    end
    
elseif plotfigs==0
    if invpar.nRuns > 1 % Added switch in case only 1 run required (ie for forward models)
        if strcmp(invpar.model{1},'FHIN') || strcmp(invpar.model{1},'XHIN')
            [Mw]=calculateMoment(invResults.mKeep);
            invResults.mKeep(end,:) = Mw;            
            model.parName{variable_param(end)} = 'Magnitude';
            
            figure('Position', [1, 1, 1200, 1000],'Visible',vis);
            
            [~,~,~,~,pax]=plotmatrix_lower(invResults.mKeep(variable_param,burning:end-blankCells)','contour',model.parName(variable_param));
            for ii=1:numel(variable_param)
                tmp_ax=get(pax(ii),'Xlabel');
                tmp_ax.String=model.parName(variable_param(ii));
            end
            
        else
            try
            figure('Position', [1, 1, 1200, 1000],'Visible',vis);
            [~,~,~,~,pax]=plotmatrix_lower(invResults.mKeep(variable_param(1:end-1),burning:end-blankCells)','contour',model.parName(variable_param(1:end-1)));
            for ii=1:numel(variable_param)-1
                tmp_ax=get(pax(ii),'Xlabel');
                tmp_ax.String=model.parName(variable_param(ii));
            end
            catch
            end
        end
        img1 = getframe(gcf);
        imwrite(img1.cdata,[outputDir,'/Figures/JointProbabilities_',num2str(burning),'.png']);
        
        % Add image to html report
        fprintf(fidHTML, '%s\r\n', '<hr>');
        fprintf(fidHTML, '%s\r\n', '<H3>Joint probabilities</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/JointProbabilities_',num2str(burning),'.png','" alt="HTML5 Icon">']);
    end
    
    %% Plot comparison betweem data, model, and residual
    
    if exist('gps','var')
        plotGPSDataModel(gps,geo,invpar, invResults, modelInput, saveName, 'y',vis)
        
        % Add image to html report
        fprintf(fidHTML, '%s\r\n', '<hr>');
        fprintf(fidHTML, '%s\r\n', '<H3>Comparison GPS Data vs. Model</H3>');
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Data_Model_horizontal.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Data_Model_vertical.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Resid_east.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Resid_north.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Model_vertical.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n', ['<img src="Figures/GPS_Contour_Resid_vertical.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n\n', ['<img src="Figures/GPS_Data_Model_scatter.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n\n', ['<img src="Figures/Vertical_Histogram.png','" alt="HTML5 Icon">']);
        fprintf(fidHTML, '%s\r\n\n', ['<img src="Figures/GPS_Data_Model_hist.png','" alt="HTML5 Icon">']);
    end
    
    % Plot InSAR data, model, residual
    if exist('insarDataCode','var')
        plotInSARDataModelResidual(insar, geo, invpar, invResults, modelInput, saveName, fidHTML, 'y',vis)
    end
    
end

fclose(fidHTML);

fprintf('Final Report Complete\n')

if plotfigs==0
    figid=findall(groot,'type','Figure');
    n_fig=size(figid,1);
    hid_fig=[];
    for ii=1:n_fig
        if strcmp(figid(ii).Visible,'off')
            hid_fig=[hid_fig,ii];
        end
    end
    close(figid(hid_fig))
    fprintf('%.0f hidden figures closed\n',size(hid_fig,2))
end

