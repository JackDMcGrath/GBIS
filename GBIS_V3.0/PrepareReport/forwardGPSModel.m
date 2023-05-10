function UGps = forwardGPSModel(xy, invpar, invResults, modelInput)

% Function to generate forward model for GPS displacements using optimal
% source parameters
%
% Usage: UGps = forwardGPSModel(xy, invpar, invResults, modelInput)
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
%
% Email: gbis.software@gmail.com
%
% Reference: 
% Bagnardi M. & Hooper A, (2018). 
% Inversion of surface deformation data for rapid estimates of source 
% parameters and uncertainties: A Bayesian approach. Geochemistry, 
% Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
% The function may include third party software.
% =========================================================================
% Last update: 8 August, 2018
%%
nu = modelInput.nu;
xy =xy';
UGps = zeros(3,length(xy(1,:)));   % Initialise matrix of modeled displacements (3 x number of observation points)

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
            mFunc{i} = invResults.model.optimal(index1:index1+7);
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
            U = hingedFaultsWithOffset(mFunc{i},xy(1:2,:),nu);
        case 'SPLT'
            mFunc{i}=invResults.model.optimal(index1:index1+11);
            U = splitFault(mFunc{i},xy(1:2,:),nu);   % Calculate 3D displacements
        case 'BSPT'
            mFunc{i}=invResults.model.optimal(index1:index1+25);     % Select source model parameters from all
            U = splithingedFault(mFunc{i},xy(1:2,:),nu);   % Calculate 3D displacements
        case 'BACK'
            mFunc{i}=invResults.model.optimal(index1:index1+10);     % Select source model parameters from all
            U = backslip(mFunc{i},xy(1:2,:),nu);   % Calculate 3D displacements
        case 'ARCT'
            mFunc{i}=[invResults.model.optimal(index1:index1+3);0]; % Select source model parameters from all; Opening set to 0;
            [xy_dist] = dist2fault_trace(invpar.fault',xy(1:2,:)); % Distances in m
            U = zeros(3,size(xy_dist,1));
            V = arctan(xy_dist,mFunc{i}(1:4));
            U(1,:) = sind(invpar.fault_strike)*V;
            U(2,:) = cosd(invpar.fault_strike)*V;
        case 'DIST'
            mFunc{i}=[invResults.model.optimal(index1:index1+4);0]; % Select source model parameters from all; Opening set to 0;
            U = zeros(3,size(invpar.zdist,1));
            V = distributed_shear(invpar.zdist,mFunc{i}(1:5));
            U(1,:) = sind(invpar.fault_strike)*V;
            U(2,:) = cosd(invpar.fault_strike)*V;
        case 'DHIN'
            mFunc{i}=[invResults.model.optimal(index1:index1+14)]; % Select source model parameters from all
            [U,~] = HingeAndDistShear(mFunc{i},xy(1:2,:),nu,invpar.zdist,invpar.fault);   % Calculate 3D displacements
    end
    UGps = UGps + U; % Calculate total displacement from sum of displacement from each source
end



