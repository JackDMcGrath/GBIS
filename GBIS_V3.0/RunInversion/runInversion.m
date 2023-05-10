function [results,invpar] = runInversion(geo, gps, insar, invpar, model, modelInput, obs, nObs)

% Function that runs the MCMC Bayesian inversion
%
% Usage: results = runInversion(geo, gps, insar, invpar, model, modelInput, obs, nObs)
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
tic; % Start timer
global outputDir  % Set global variables

% Calculate remaining time using time since last n printouts
backTime=4;
preiKeep=zeros(1,1+backTime);
pretoc=zeros(1,1+backTime);
%% Set starting model parameters
nu = modelInput.nu; % Poisson's ratio

model.trial = model.m; % First trial using starting parameters

model.range = model.upper - model.lower; % Calculate range of m parameters

% Check that starting model is within bounds ( but only if not doing a
% forward model)
if invpar.nRuns > 1
    invparam=find(model.step);
    if sum(model.m(invparam) > model.upper(invparam))>0 || sum(model.m(invparam) < model.lower(invparam))>0
        fprintf('#\tParameter\tLower_Bound\tStart_Model\tUpper_Bound\n')
        ix = find(model.m > model.upper | model.m < model.lower);
        for ii=1:size(ix,1)
            fprintf('%d\t%s\t%.6f\t%.6f\t%.6f\n', ix(ii), model.parName{ix(ii)}, model.lower(ix(ii)), model.m(ix(ii)), model.upper(ix(ii)))
        end
        error('Starting model is out of bounds')
    end
end

nModel = length(model.m);    % Number of model parameters to invert for
probTarget = 0.5^(1/nModel); % Target probability used in sensitivity test
probSens = zeros(nModel,1);  % Initialise vector for sensitivity test

iKeep = 0; % Initialiase kept iterations counter
iReject = 0;    % Initialise rejected iterations counter
iKeepSave = iKeep;  % Initialise saving schedule for kept iterations
iRejectSave = iReject; % Initialise saving schedule for rejected iterations

mKeep= zeros(nModel,invpar.nSave,'single');   % Initialise matrix of model parameters to keep
PKeep= zeros(1,invpar.nSave,'single');   % Initialise probability vector
RKeep= zeros(1,invpar.nSave,'single');   % Initialise RMS vector

POpt = -1e99; % Set initial optimal probability
T = invpar.TSchedule(1); % Set first temperature
iTemp = 0;  % Initialise temperature schedule index (for initial Simulated Annealing)
nTemp = length(invpar.TSchedule); % Number of temperature steps (for initial Simulated Annealing)

%% Start core of inversion

sensitivityTest = 0; % Switch off sensitivity test at first iteration

while iKeep < invpar.nRuns  % While number of iterations is < than number of runs...
    if iKeep/invpar.TRuns == round(iKeep/invpar.TRuns) & iTemp < nTemp % Follow temperature schedule
        iTemp = iTemp + 1;
        T = invpar.TSchedule(iTemp);    % Assign temperature from T schedule
        if iKeep > 0
            model.trial = model.optimal;
        end
        if T ==1                        % Set Hyperparameter when T reaches 1 (Hyperparameter is currently not is use!!!)
            setHyperParameter = 1;
        else
            setHyperParameter = 0;
        end
    end
    
    
    if sum(iKeep == invpar.sensitivitySchedule)>0   % Check if it's time for sensitivity test based on schedule
        sensitivityTest = 1; % Switch on sensitivity test
    end
    
    %% Calculate 3D displacements from model
    
    UTot = zeros(3,nObs);   % Initialise matrix of modeled displacements (3 x number of observation points)
    [modelTypes,ix,IX]=unique(invpar.model);
    modelTypes=modelTypes(IX(ix)); % Keep models in correct order
    mCount=1; % Number of source models checked
    for i = 1:size(modelTypes,2) % For each type of source model...
        sourceModel=modelTypes{i};
        modelIx=find(strcmp(invpar.model,modelTypes{i}));
        for ii = 1:length(modelIx);
            ix=modelIx(ii);
            index1 = model.mIx(ix);
            switch invpar.model{ix}
                case 'MOGI'
                    mFunc{mCount} = model.trial(index1:index1+3);    % Select source model parameters from all
                    U = mogi(mFunc{mCount},obs,nu);               % Calculate 3D displacements
                case 'MCTG'
                    mFunc{mCount} = model.trial(index1:index1+4);    % Select source model parameters from all
                    U = mctigueSource(mFunc{mCount},obs(1:2,:),nu);        % Calculate 3D displacements
                case 'YANG'
                    mFunc{mCount} = model.trial(index1:index1+7);    % Select source model parameters from all
                    U = yangSource(mFunc{mCount},obs,nu);            % Calculate 3D displacements
                case 'PENN'
                    mFunc{mCount}=model.trial(index1:index1+4);      % Select source model parameters from all
                    U = pennySource(mFunc{mCount},obs,nu);           % Calculate 3D displacements
                case 'SILL'
                    mFunc{mCount}=[model.trial(index1:index1+2); 0; model.trial(index1+3:index1+5); 0; 0; model.trial(index1+6)]; % Select source model parameters from all; Dip set to 0; Slip set to 0;
                    U = disloc(mFunc{mCount},obs(1:2,:),nu);
                case 'DIKE'
                    mFunc{mCount}=[model.trial(index1:index1+6); 0; 0; model.trial(index1+7)]; % Select source model parameters from all; Slip set to 0;
                    U = disloc(mFunc{mCount},obs(1:2,:),nu);
                case 'FAUL'
                    mFunc{mCount}=[model.trial(index1:index1+8);0]; % Select source model parameters from all; Opening set to 0;
                    %                 trace.linktotrace=modelInput.fault{i}.linktotrace; % Flag to ensure location of faultsplit is kept onto surface trace if faulttracefile is provided
                    %                 if trace.linktotrace==1 && isfield(geo,'faulttracefile') % If flagged, ensure the split is kept to the surface trace
                    %                     trace.tracefile=load(geo.faulttracefile);
                    %                     trace.tracefile = llh2local([trace.tracefile'; zeros(1,size(trace.tracefile,1))], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
                    %                     try
                    %                         [mFunc{i}(6),mFunc{i}(7)]=point_to_lineXY(mFunc{i}(6),mFunc{i}(7),trace.tracefile);
                    %                         model.trial(index1+5:index1+6)=mFunc{i}(6:7);
                    %                     catch
                    %                         fprintf('Error locating split on trace for run %.0f\n',iKeep)
                    %                     end
                    %                 elseif trace.linktotrace==1
                    %                     trace.linktotrace=0;
                    %                     fprintf('No Fault Trace File Provided. Split location unconstrained for SPLIT %.0f',i)
                    %                 end
                    U = disloc(mFunc{mCount},obs(1:2,:),nu);
                case 'ARCT'
                    mFunc{mCount}=[model.trial(index1:index1+3);0]; % Select source model parameters from all; Opening set to 0;
                    U = zeros(3,nObs);
                    V = arctan(invpar.dist,mFunc{mCount}(1:4));
                    U(1,:) = sind(geo.fault_strike)*V;
                    U(2,:) = cosd(geo.fault_strike)*V;
                case 'DIST'
                    mFunc{mCount}=[model.trial(index1:index1+4);0]; % Select source model parameters from all; Opening set to 0;
                    U = zeros(3,nObs);
                    V = distributed_shear(invpar.zdist,mFunc{mCount}(1:5));
                    U(1,:) = sind(geo.fault_strike)*V;
                    U(2,:) = cosd(geo.fault_strike)*V;
                case 'BACK'
                    mFunc{mCount}=[model.trial(index1:index1+10);0]; % Select source model parameters from all; Opening set to 0;
                    U = backslip(mFunc{mCount},obs(1:2,:),nu);
                case 'HING'
                    mFunc{mCount}=model.trial(index1:index1+10);     % Select source model parameters from all
                    U = hingedDikes(mFunc{mCount},obs(1:2,:),nu);    % Calculate 3D displacements
                case 'FHIN'
                    model.trial(index1:index1+13)=linked_paramaters(model.trial(index1:index1+13),modelInput.fhin{i}.link);
                    mFunc{mCount}=model.trial(index1:index1+13);     % Select source model parameters from all
                    U = hingedFaults(mFunc{mCount},obs(1:2,:),nu);   % Calculate 3D displacements
                case 'XHIN'
                    model.trial(index1:index1+14)=linked_paramaters(model.trial(index1:index1+14),modelInput.xhin{ii}.link);
                    mFunc{mCount}=model.trial(index1:index1+14);     % Select source model parameters from all
                    U = hingedFaultsWithOffset(mFunc{mCount},obs(1:2,:),nu);   % Calculate 3D displacements
                case 'SPLT'
                    mFunc{mCount}=model.trial(index1:index1+11);     % Select source model parameters from all
                    trace.linktotrace=modelInput.split{ix}.linktotrace; % Flag to ensure location of faultsplit is kept onto surface trace if faulttracefile is provided
                    if trace.linktotrace==1 && isfield(geo,'faulttracefile') % If flagged, ensure the split is kept to the surface trace
                        trace.tracefile=load(geo.faulttracefile);
                        trace.tracefile = llh2local([trace.tracefile'; zeros(1,size(trace.tracefile,1))], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
                        try
                            [mFunc{mCount}(6),mFunc{mCount}(7)]=point_to_lineXY(mFunc{mCount}(6),mFunc{mCount}(7),trace.tracefile);
                            model.trial(index1+5:index1+6)=mFunc{mCount}(6:7);
                        catch
                            fprintf('Error locating split on trace for run %.0f\n',iKeep)
                        end
                    elseif trace.linktotrace==1
                        trace.linktotrace=0;
                        fprintf('No Fault Trace File Provided. Split location unconstrained for SPLIT %.0f',i)
                    end
                    U = splitFault(mFunc{mCount},obs(1:2,:),nu,trace);   % Calculate 3D displacements
                case 'BSPT'
                    model.trial(index1:index1+25)=linked_paramaters(model.trial(index1:index1+25),modelInput.bsplt{ix}.link);
                    mFunc{mCount}=model.trial(index1:index1+25);     % Select source model parameters from all
                    trace.linktotrace=modelInput.bsplt{ix}.linktotrace; % Flag to ensure location of faultsplit is kept onto surface trace if faulttracefile is provided
                    if trace.linktotrace==1 && isfield(geo,'faulttracefile') % If flagged, ensure the split is kept to the surface trace
                        trace.tracefile=load(geo.faulttracefile);
                        trace.tracefile = llh2local([trace.tracefile'; zeros(1,size(trace.tracefile,1))], geo.referencePoint)*1000; % Convert geographic coordinates to local cooridinates
                        try
                            [mFunc{mCount}(2),mFunc{mCount}(3)]=point_to_lineXY(mFunc{mCount}(2),mFunc{mCount}(3),trace.tracefile);
                            model.trial(index1+1:index1+2)=mFunc{mCount}(2:3);
                        catch
                            fprintf('Error locating split on trace for run %.0f\n',iKeep)
                        end
                    elseif trace.linktotrace==1
                        trace.linktotrace=0;
                        fprintf('No Fault Trace File Provided. Split location unconstrained for BSPLT %.0f',i)
                    end
                    U = splithingedFault(mFunc{mCount},obs(1:2,:),nu);   % Calculate 3D displacements
                case 'DHIN'
                    mFunc{mCount}=model.trial(index1:index1+14); % Select source model parameters from all
                    [U, right] = HingeAndDistShear(mFunc{mCount},obs(1:2,:),nu,invpar.zdist,invpar.fault);   % Calculate 3D displacements
                    invpar.right=right;
            end
            UTot = UTot + U; % Calculate total displacement from sum of displacement from each source
            mCount = mCount+1;
        end
    end
    
    insarParIx = model.mIx(end); % Identify first model parameter not related to source model (e.g., InSAR offset, ramp, etc.)
    
    U = UTot; % Reassign U to Utot for simplicity
    
    %% Convert 3D displacement to LOS displacement and calculate residuals
    resExp = 0; % Initialise (Gm - d) * InvCov * (Gm - d)'
    
    if ~isempty(insar)
        
        ULos = []; % Initialise line-of-sight Gm vector
        
        for j = 1 : length(insar)
            UEast = -cosd(insar{j}.dHeading).* sind(insar{j}.dIncidence); % East unit vector
            UNorth = sind(insar{j}.dHeading).* sind(insar{j}.dIncidence); % North unit vector
            UVert = cosd(insar{j}.dIncidence); % Vertical unit vector
            
            ULos{j} = UEast.* U(1,insar{j}.ix) + ...
                UNorth.* U(2,insar{j}.ix) + ...             % Convert to line of sight displacement
                UVert.* U(3,insar{j}.ix);
            
            if insar{j}.constOffset == 'y'
                ULos{j} = ULos{j} + model.trial(insarParIx);  % Add constant offset
                
                insarParIx = insarParIx + 1; % Change model parameter index for next step
            end
            
            if insar{j}.rampFlag == 'y'
                ULos{j} = ULos{j} + model.trial(insarParIx)*insar{j}.obs(:,1)' + ...
                    model.trial(insarParIx+1)*insar{j}.obs(:,2)'; % Add linear ramp
                
                insarParIx = insarParIx + 2; % Change model parameter index for next step if necessary
            end
            
            resInsar{j} = (ULos{j} - insar{j}.dLos); % Calculate (Gm - d), residuals
            resExp = resExp + resInsar{j}* insar{j}.invCov* resInsar{j}'; % (Gm - d) * InvCov * (Gm - d)'
        end
    end
    
    
    %% Calculate GPS residuals
    
    if ~isempty(gps)
        for i = 1 : length(gps)
            rGps = (U(1:3,gps{i}.ix) - gps{i}.displacements(1:3,:));%.*gps.enu_weight; % Residual GPS displacement
            resExp = resExp + rGps(:)' * gps{i}.invCov * rGps(:) * gps{i}.weight; % (Gm - d) * InvCov * (Gm - d)'
        end
    end
    
    %% Continue inversion ...
    
    if setHyperParameter == 1
        %hyperPrev = resExp/nObs; % set hyperparameter on first reaching T=1;
        hyperPrev = 1; % set hyperparameter to 1;
        model.trial(end) = log10(hyperPrev);
        setHyperParameter = 0; % Switch setHyperParameter off
    end
    
    if isempty(insar)
        hyperParam = 1;
    else
        hyperParam = 1;
        %hyperParam = 10^model.trial(end);
    end
    
    % !! Currently hyperparameter is set to 1
    P = -resExp/(2*hyperParam); % Probability is exp of P
    
    if iKeep>0
        PRatio = (hyperPrev/hyperParam)^(nObs/2)*exp((P-PPrev)/T);  % Probability ratio
    else
        PRatio=1; % Set to 1 for first iteration (always keep first iteration)
    end
    
    %% Perform sensitivity test if necessary and change step size
    
    if sensitivityTest > 1
        probSens(sensitivityTest-1) = PRatio; % Assign probability to current model parameter
        if sensitivityTest > nModel % Check if sensitivity test has finished
            if iKeepSave > 0
                rejectionRatio = (iReject - iRejectSave)/(iKeep - iKeepSave); % Calculate rejection rate
                probTarget = probTarget * rejectionRatio * 1/0.77; % Adjust target probability to reach 77% rejection rate
                probTarget(probTarget<1e-06) = 1e-06; % Prevent from reaching zero.
            end
            sensitivityTest = 0;    % Swtich off sensitivity test
            probSens(probSens > 1) = 1./probSens(probSens > 1);
            PDiff = probTarget - probSens;
            indexP = PDiff > 0; % Select model parameters for which to adjust model step
            model.step(indexP) = model.step(indexP).*exp(-PDiff(indexP)/probTarget*2);  % assign new model step
            indexP = PDiff < 0; % Select remaining model parameters
            model.step(indexP) = model.step(indexP).*exp(-PDiff(indexP)/(1-probTarget)*2); % assign new model step
            model.step(model.step > model.range) = model.range(model.step > model.range); % Check if step is within range
            iKeepSave = iKeep;
            iRejectSave = iReject;
        end
        
    else
        iKeep = iKeep + 1;
        if PRatio >= rand(1,1)  % If condions are met, keep model trial
            model.m = model.trial; % Substitute m with model trial
            mKeep(:,iKeep) = model.m;   % Keep model trial
            PKeep(:,iKeep) = P;         % P of current model
            RKeep(:,iKeep) = resExp;   % resExp of current model
            
            PPrev = P;  % Assign current P to PPrev for next trial
            hyperPrev = hyperParam; % Assign current Hyperparameter for next trial
            
            if -resExp > POpt   % Update current optimal model if likelihood is higher
                model.optimal = model.m;
                model.funcOpt = mFunc;
                POpt = -resExp;
            end
        else                    % Reject model trial and keep previous model
            iReject = iReject + 1;
            mKeep(:,iKeep) = mKeep(:,iKeep-1);
            PKeep(:,iKeep) = PKeep(:,iKeep-1);
            RKeep(:,iKeep) = RKeep(:,iKeep-1);
        end
        
        if iKeep/invpar.nSave == round(iKeep/invpar.nSave) % display and save results at regular intervals (1000 or 10000 iterations)
            if iKeep >= 20000           % Increase time step for saving/displaying after 20000 iterations
                invpar.nSave = 10000;
            end
            
            % Print current status of inversion to screen
            disp('=========================================================')
            disp(['Model: ',invpar.model{:}, ' (',outputDir, ')'])
            disp([num2str(iKeep),'/',num2str(invpar.nRuns),' model trials in ', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'),' (ETC: ',datestr(datenum(0,0,0,0,0,toc+((toc-pretoc(end-backTime))*(invpar.nRuns-iKeep)/(iKeep-preiKeep(end-backTime)))),'HH:MM:SS'),')'])
            preiKeep=[preiKeep,iKeep];
            pretoc=[pretoc,toc];
            disp(['Optimal Prob = exp(',num2str(POpt),')'])
            disp(['Hyperparameter=',num2str(hyperParam)])
            disp([num2str(iReject),' models rejected:', num2str((iReject/iKeep)*100),'% of model trials.'])
            
            % allocate space for next blocks to keep
            mKeep(:,iKeep + invpar.nSave) = 0;
            PKeep(:,iKeep + invpar.nSave) = 0;
            RKeep(:,iKeep + invpar.nSave) = 0;
            
            % Save results to temporary file for insepction during
            % inversion
            save([outputDir,'/temporary.mat'], 'geo', 'mKeep', 'PKeep', 'model', 'gps', 'insar', 'invpar', 'geo', 'modelInput','-v7.3');
            
            
            % Display current optimal model parameters on screen
            for i=1:length(invpar.model)
                if invpar.model{i} == 'MOGI'
                    fprintf('MOGI center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('MOGI center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('MOGI depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('MOGI volume change: %f\n',(model.funcOpt{i}(4,:)));
                elseif invpar.model{i} == 'YANG'
                    fprintf('YANG centroid X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('YANG centroid Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('YANG centroid depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('YANG major axis: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('YANG minor axis: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('YANG majax strike: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('YANG majax plunge: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('YANG DP/mu: %f\n',(model.funcOpt{i}(8,:)));
                elseif invpar.model{i} == 'MCTG'
                    fprintf('MCTG center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('MCTG center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('MCTG center depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('MCTG radius: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('MCTG DP/mu: %f\n',(model.funcOpt{i}(5,:)));
                elseif invpar.model{i} == 'PENN'
                    fprintf('PENNY center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('PENNY center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('PENNY center depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('PENNY radius: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('PENNY DP/mu: %f\n',(model.funcOpt{i}(5,:)));
                elseif invpar.model{i} == 'SILL'
                    fprintf('SILL length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('SILL width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('SILL depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('SILL strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('SILL edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('SILL edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('SILL opening: %f\n',(model.funcOpt{i}(10,:)));
                elseif invpar.model{i} == 'DIKE'
                    fprintf('DIKE length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('DIKE width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('DIKE depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('DIKE dip: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('DIKE strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('DIKE edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('DIKE edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('DIKE opening: %f\n',(model.funcOpt{i}(10,:)));
                elseif invpar.model{i} == 'FAUL'
                    fprintf('FAULT length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('FAULT width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('FAULT depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('FAULT dip: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('FAULT strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('FAULT edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('FAULT edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('FAULT strike-slip component: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('FAULT dip-slip component: %f\n',(model.funcOpt{i}(9,:)));
                elseif invpar.model{i} == 'HING'
                    fprintf('DIKE1 edge center X: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('DIKE1 edge center Y: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('DIKE1 length: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('DIKE1 width: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('DIKE1 depth: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('DIKE1 dip: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('DIKE1 opening: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('DIKE2 width: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('DIKE2 dip: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('DIKE2 opening: %f\n',(model.funcOpt{i}(10,:)));
                    fprintf('Strike: %f\n',(model.funcOpt{i}(11,:)));
                elseif invpar.model{i} == 'FHIN'
                    fprintf('FAULT length: %.2f km\n',(model.funcOpt{i}(1,:)*1e-3));
                    fprintf('FAULT1 width: %.2f km\n',(model.funcOpt{i}(2,:)*1e-3));
                    fprintf('FAULT1 depth: %.2f km\n',(model.funcOpt{i}(3,:)*1e-3));
                    fprintf('FAULT1 dip: %.1f deg\n',(model.funcOpt{i}(4,:)));
                    fprintf('FAULT strike: %.0f deg\n',(model.funcOpt{i}(5,:)));
                    fprintf('FAULT edge center X: %.0f\n',(model.funcOpt{i}(6,:)));
                    fprintf('FAULT edge center Y: %.0f\n',(model.funcOpt{i}(7,:)));
                    fprintf('FAULT1 strike-slip component: %.2f mm\n',(model.funcOpt{i}(8,:)*1e3));
                    fprintf('FAULT1 dip-slip component: %.2f mm\n',(model.funcOpt{i}(9,:)*1e3));
                    fprintf('FAULT2 width: %.2f km\n',(model.funcOpt{i}(10,:)*1e-3));
                    fprintf('FAULT2 dip: %.1f deg\n',(model.funcOpt{i}(11,:)));
                    fprintf('FAULT2 strike-slip component: %.2f mm\n',(model.funcOpt{i}(12,:)*1e3));
                    fprintf('FAULT2 dip-slip component: %.2f mm\n',(model.funcOpt{i}(13,:)*1e3));
                elseif invpar.model{i} == 'XHIN'
                    fprintf('FAULT length: %.2f km\n',(model.funcOpt{i}(1,:)*1e-3));
                    fprintf('FAULT1 width: %.2f km\n',(model.funcOpt{i}(2,:)*1e-3));
                    fprintf('FAULT1 depth: %.2f km\n',(model.funcOpt{i}(3,:)*1e-3));
                    fprintf('FAULT1 dip: %.1f deg\n',(model.funcOpt{i}(4,:)));
                    fprintf('FAULT strike: %.0f deg\n',(model.funcOpt{i}(5,:)));
                    fprintf('FAULT edge center X: %.0f\n',(model.funcOpt{i}(6,:)));
                    fprintf('FAULT edge center Y: %.0f\n',(model.funcOpt{i}(7,:)));
                    fprintf('FAULT1 strike-slip component: %.2f mm\n',(model.funcOpt{i}(8,:)*1e3));
                    fprintf('FAULT1 dip-slip component: %.2f mm\n',(model.funcOpt{i}(9,:)*1e3));
                    fprintf('FAULT2 width: %.2f km\n',(model.funcOpt{i}(10,:)*1e-3));
                    fprintf('FAULT2 dip: %.1f deg\n',(model.funcOpt{i}(11,:)));
                    fprintf('FAULT2 strike-slip component: %.2f mm\n',(model.funcOpt{i}(12,:)*1e3));
                    fprintf('FAULT2 dip-slip component: %.2f mm\n',(model.funcOpt{i}(13,:)*1e3));
                    fprintf('FAULT Perpendicular offset: %.2f km\n',(model.funcOpt{i}(14,:)*1e-3));
                elseif invpar.model{i} == 'SPLT'
                    fprintf('SPLIT length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('SPLIT base: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('SPLIT locking depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('SPLIT strikesuch that the mean chance of acce: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('SPLIT X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('SPLIT Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('SPLIT1 dip: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('SPLIT1 strike-slip component: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('SPLIT1 dip-slip component: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('SPLIT2 dip: %f\n',(model.funcOpt{i}(10,:)));
                    fprintf('SPLIT2 strike-slip component: %f\n',(model.funcOpt{i}(11,:)));
                    fprintf('SPLIT2 dip-slip component: %f\n',(model.funcOpt{i}(12,:)));
                elseif invpar.model{i} == 'BSPT'
                    fprintf('BSPLT Strike: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('BSPLT X: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('BSPLT Y: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('BSPLT1 length: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('BSPLT1 locking depth: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('BSPLT1 base: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('BFAULT1 width: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('BFAULT1 dip: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('BFAULT1 strike-slip component: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('BFAULT1 dip-such that the mean chance of acceslip component: %f\n',(model.funcOpt{i}(10,:)));
                    fprintf('BFAULT2 width: %f\n',(model.funcOpt{i}(11,:)));
                    fprintf('BFAULT2 dip: %f\n',(model.funcOpt{i}(12,:)));
                    fprintf('BFAULT2 strike-slip component: %f\n',(model.funcOpt{i}(13,:)));
                    fprintf('BFAULT2 dip-slip component: %f\n',(model.funcOpt{i}(14,:)));
                    fprintf('BSPLT2 length: %f\n',(model.funcOpt{i}(15,:)));
                    fprintf('BSPLT2 locking depth: %f\n',(model.funcOpt{i}(16,:)));
                    fprintf('BSPLT2 base: %f\n',(model.funcOpt{i}(17,:)));
                    fprintf('BFAULT3 width: %f\n',(model.funcOpt{i}(18,:)));
                    fprintf('BFAULT3 dip: %f\n',(model.funcOpt{i}(19,:)));
                    fprintf('BFAULT3 strike-slip component: %f\n',(model.funcOpt{i}(20,:)));
                    fprintf('BFAULT3 dip-slip component: %f\n',(model.funcOpt{i}(21,:)));
                    fprintf('BFAULT4 width: %f\n',(model.funcOpt{i}(22,:)));
                    fprintf('BFAULT4 dip: %f\n',(model.funcOpt{i}(23,:)));
                    fprintf('BFAULT4 strike-slip component: %f\n',(model.funcOpt{i}(24,:)));
                    fprintf('BFAULT4 dip-slip component: %f\n',(model.funcOpt{i}(25,:)));
                elseif invpar.model{i} == 'BACK'
                    fprintf('BACKSLIP length: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('BACKSLIP width: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('BACKSLIP depth: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('BACKSLIP dip: %f\n',(model.funcOpt{i}(4,:)));
                    fprintf('BACKSLIP strike: %f\n',(model.funcOpt{i}(5,:)));
                    fprintf('BACKSLIP edge center X: %f\n',(model.funcOpt{i}(6,:)));
                    fprintf('BACKSLIP edge center Y: %f\n',(model.funcOpt{i}(7,:)));
                    fprintf('BACKSLIP strike-slip component: %f\n',(model.funcOpt{i}(8,:)));
                    fprintf('BACKSLIP dip-slip component: %f\n',(model.funcOpt{i}(9,:)));
                    fprintf('BACKSLIP Lvel: %f\n',(model.funcOpt{i}(10,:)));
                    fprintf('BACKSLIP Rvel: %f\n',(model.funcOpt{i}(11,:)));
                elseif invpar.model{i} == 'ARCT'
                    fprintf('ARCTAN Strike-Slip: %f\n',(model.funcOpt{i}(1,:)));
                    fprintf('ARCTAN Locking Depth: %f\n',(model.funcOpt{i}(2,:)));
                    fprintf('ARCTAN Assymmetry: %f\n',(model.funcOpt{i}(3,:)));
                    fprintf('ARCTAN Horizontal Offset: %f\n',(model.funcOpt{i}(4,:)));
                elseif invpar.model{i} == 'DIST'
                    fprintf('DISTRIBUTED SHEAR Strike-Slip: %.2f mm/yr\n',(model.funcOpt{i}(1,:)));
                    fprintf('DISTRIBUTED SHEAR Locking Depth: %.2f km\n',(model.funcOpt{i}(2,:)));
                    fprintf('DISTRIBUTED SHEAR Center: %.2f km\n',(model.funcOpt{i}(3,:)));
                    fprintf('DISTRIBUTED SHEAR Width: %.2f km\n',(model.funcOpt{i}(4,:)));
                elseif invpar.model{i} == 'DHIN'
                    fprintf('FAULT length: %.2f km\n',(model.funcOpt{i}(1,:)*1e-3));
                    fprintf('FAULT1 width: %.2f km\n',(model.funcOpt{i}(2,:)*1e-3));
                    fprintf('FAULT1 depth: %.2f km\n',(model.funcOpt{i}(3,:)*1e-3));
                    fprintf('FAULT1 dip: %.1f deg\n',(model.funcOpt{i}(4,:)));
                    fprintf('FAULT strike: %.0f deg\n',(model.funcOpt{i}(5,:)));
                    fprintf('FAULT edge center X: %.0f\n',(model.funcOpt{i}(6,:)));
                    fprintf('FAULT edge center Y: %.0f\n',(model.funcOpt{i}(7,:)));
                    fprintf('FAULT1 strike-slip component: %.2f mm\n',(model.funcOpt{i}(8,:)*1e3));
                    fprintf('FAULT1 dip-slip component: %.2f mm\n',(model.funcOpt{i}(9,:)*1e3));
                    fprintf('FAULT2 width: %.2f km\n',(model.funcOpt{i}(10,:)*1e-3));
                    fprintf('FAULT2 dip: %.1f deg\n',(model.funcOpt{i}(11,:)));
                    fprintf('FAULT2 strike-slip component: %.2f mm\n',(model.funcOpt{i}(12,:)*1e3));
                    fprintf('FAULT2 dip-slip component: %.2f mm\n',(model.funcOpt{i}(13,:)*1e3));
                    fprintf('DISTRIBUTED SHEAR Strike-Slip: %.2f mm/yr\n',(model.funcOpt{i}(14,:)));
                    fprintf('DISTRIBUTED SHEAR Depth: %.2f km\n',(model.funcOpt{i}(3,:)*1e-3 + sind(-model.funcOpt{i}(4,:))*model.funcOpt{i}(2,:)*1e-3));
                    fprintf('DISTRIBUTED SHEAR Left Edge: %.2f km\n',(model.funcOpt{i}(15,:)));
                    fprintf('DISTRIBUTED SHEAR Right Edge: %.2f km\n',invpar.right);
                end
            end
        end
    end
    
    
    
    if sensitivityTest > 0  % Perform sensitivity test (no models are kept during this phase!)
        randomStep = zeros(nModel,1);
        randomStep(sensitivityTest) = model.step(sensitivityTest) * sign(randn(1,1))/2; % Assign random step
        model.trial = model.m + randomStep; % New model trial
        % Check that new model trial is withing bounds
        if model.trial(sensitivityTest) > model.upper(sensitivityTest)
            model.trial(sensitivityTest) = model.trial(sensitivityTest) - model.step(sensitivityTest);
        end
        
        hyperParam = hyperPrev;
        sensitivityTest = sensitivityTest + 1; % Move index to that of next parameter until all parameters are done
    else
        randomStep = model.step.*(rand(nModel,1)-0.5)*2;     % Make random step
        model.trial = model.m + randomStep;                 % Assign new model trial to previous + random step
        % Check that new model trial is withing bounds
        model.trial(model.trial > model.upper) = 2 * model.upper(model.trial > model.upper) - ...
            model.trial(model.trial > model.upper);
        
        model.trial(model.trial < model.lower) = 2 * model.lower(model.trial < model.lower) - ...
            model.trial(model.trial < model.lower);
    end
    
    %     if iKeep > (0.1 * invpar.nRuns)
    %         test_param = find(model.step ~= 0);
    %         n_test = numel(test_param);
    %         search_range = 0.05 * (model.upper(test_param) - model.lower(test_param));
    %         max_lim = model.upper(test_param)-search_range;
    %         min_lim = model.lower(test_param)+search_range;
    %         model_range = range(mKeep(test_param,iKeep-0.1*invpar.nRuns:iKeep)')';
    %         converged = test_param(find(model_range<search_range));
    %         if ~isempty(converged)
    %             mean_conv = mean(mKeep(converged,iKeep-0.1*invpar.nRuns:iKeep)')';
    %         end
    %     end
end


%% Clean up and prepare results
if invpar.nRuns > 1 % Added switch in case only 1 run required (ie for forward models)
    mKeep(:, end - invpar.nSave) = []; % Remove unused preallocated memory
    PKeep(:, end - invpar.nSave) = []; % Remove unused preallocated memory
    RKeep(:, end - invpar.nSave) = []; % Remove unused preallocated memory
end

results.mKeep = mKeep;
results.PKeep = PKeep;
results.RKeep = RKeep;
results.model = model;
results.optimalmodel = model.funcOpt;
results.reject_percent = (iReject/iKeep)*100;

