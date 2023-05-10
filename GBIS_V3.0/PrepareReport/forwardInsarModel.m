function ULos = forwardInSARModel(insar,xy,invpar,invResults,modelInput,geo,Heading,Inc,constOffset,xRamp,yRamp)

% Function to generate forward model for InSAR displacements using optimal
% source parameters
%
% Usage: ULos = forwardInSARModel(insar,xy,invpar,invResults,modelInput,geo,Heading,Inc,constOffset,xRamp,yRamp)
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
nu = modelInput.nu;
xy = [xy(:,2:3)'; zeros(1,length(xy(:,1)))];
UTot = zeros(3,length(xy(1,:)));   % Initialise matrix of modeled displacements (3 x number of observation points)

for i = 1:invpar.nModels % For each source model...
    index1 = invResults.model.mIx(i);
    switch invpar.model{i}
        case 'MOGI'
            mFunc{i} = invResults.model.optimal(index1:index1+3); % Parameters to invert for.
            U = mogi(mFunc{i},xy,nu);
        case 'MCTG'
            mFunc{i} = invResults.model.optimal(index1:index1+4);
            U = mctigueSource(mFunc{i},xy(1:2,:),nu);
        case 'YANG'
            mFunc{i} = invResults.model.optimal(index1:iindex1+7);
            U = yangSource(mFunc{i},xy,nu);
        case 'PENN'
            mFunc{i}= invResults.model.optimal(index1:index1+4);
            U = pennySource(mFunc{i},xy,nu);
        case 'SILL'
            mFunc{i}=[invResults.model.optimal(index1:index1+2);0; ...
                invResults.model.optimal(index1+3:index1+5);0;0; ...
                invResults.model.optimal(index1+6)]; 
            U = disloc(mFunc{i},xy(1:2,:),nu); 
        case 'DIKE'
            mFunc{i}=[invResults.model.optimal(index1:index1+6);0;0;...
                invResults.model.optimal(index1+7)]; 
            U = disloc(mFunc{i},xy(1:2,:),nu); 
        case 'FAUL'
            mFunc{i}=[invResults.model.optimal(index1:index1+8);0]; 
            U = disloc(mFunc{i},xy(1:2,:),nu); 
        case 'HING'
            mFunc{i}=invResults.model.optimal(index1:index1+10);
            U = hingedDikes(mFunc{i},xy(1:2,:),nu);
        case 'FHIN'
            mFunc{i}=invResults.model.optimal(index1:index1+12);
            U = hingedFaults(mFunc{i},xy(1:2,:),nu);
        case 'XHIN'
            mFunc{i}=invResults.model.optimal(index1:index1+13);
            U = hingedFaults(mFunc{i},xy(1:2,:),nu);
        case 'ARCT'
            mFunc{i}=[invResults.model.optimal(index1:index1+3);0]; % Select source model parameters from all; Opening set to 0;
            [invpar.dist] = dist2fault_trace(invpar.fault',xy([1 2],:)); % Distances in m
            U = zeros(3,size(invpar.dist,1));
            V = arctan(invpar.dist,mFunc{i}(1:4));
            U(1,:) = sind(invpar.fault_strike)*V;
            U(2,:) = cosd(invpar.fault_strike)*V;
    end
    UTot = UTot + U; % Calculate total displacement from sum of displacement from each source
end
insarParIx = invResults.model.mIx(end); % Identify first model parameter not related to source model (e.g., offset, ramp, etc.)

ULos = [];

UEast = -cosd(Heading).* sind(Inc); % East unit vector
UNorth = sind(Heading).* sind(Inc); % North unit vector
UVert = cosd(Inc); % Vertical unit vector

ULos = UEast'.* UTot(1,:) + ...
    UNorth'.* UTot(2,:) + ... % Convert to line of sight displacement
    UVert'.* UTot(3,:);

if insar.constOffset == 'y'
    ULos = ULos + invResults.model.optimal(constOffset);  % Add constant offset
end

if insar.rampFlag == 'y'
    ULos = ULos + invResults.model.optimal(xRamp)*xy(1,:) + ...
        invResults.model.optimal(yRamp)*xy(2,:); % Add ramp
end

end
