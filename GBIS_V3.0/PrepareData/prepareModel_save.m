function model = prepareModel(modelInput, invpar, insar, gps)

% Function to prepare model parameters
%
% Usage: model = prepareModel(modelInput, invpar, insar, gps)
% Input Parameters:
%       modelInput: parameters read from input file
%       invpar: inversion parameters
%
% Output Parameters:
%       model: structure containing model settings to be used for inversion
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
% Last update: 02/05/2017

%% Initialize model vectors
mIx = zeros(invpar.nModels+1, 1);
mIx(1) = 1;
model.m = zeros(500,1);
model.step = model.m;
model.lower = model.m;
model.upper = model.m;

%% Assign model parameters from input file
for i = 1 : invpar.nModels
    index1 = mIx(i);
    switch invpar.model{i}
        case 'MOGI'
            nParameters = 4;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'MOGI'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.mogi{ii}.start;
            model.step(index1:index2) = modelInput.mogi{ii}.step;
            model.lower(index1:index2) = modelInput.mogi{ii}.lower;
            model.upper(index1:index2) = modelInput.mogi{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'MOGI X'; 'MOGI Y'; 'MOGI Depth'; 'MOGI DV'};
        case 'YANG'
            nParameters = 8;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'YANG'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.yang{ii}.start;
            model.step(index1:index2) = modelInput.yang{ii}.step;
            model.lower(index1:index2) = modelInput.yang{ii}.lower;
            model.upper(index1:index2) = modelInput.yang{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'YANG X'; 'YANG Y'; 'YANG Depth'; 'YANG majAx'; 'YANG a/r'; ...
                'YANG strike'; 'YANG Plunge'; 'YANG DP/mu'};
        case 'MCTG'
            nParameters = 5;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'MCTG'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.mctigue{ii}.start;
            model.step(index1:index2) = modelInput.mctigue{ii}.step;
            model.lower(index1:index2) = modelInput.mctigue{ii}.lower;
            model.upper(index1:index2) = modelInput.mctigue{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'MCTG X'; 'MCTG Y'; 'MCTG Depth'; 'MCTG Radius'; 'MCTG DP/mu'};
        case 'PENN'
            nParameters = 5;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'PENN'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.penny{ii}.start;
            model.step(index1:index2) = modelInput.penny{ii}.step;
            model.lower(index1:index2) = modelInput.penny{ii}.lower;
            model.upper(index1:index2) = modelInput.penny{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'PENN X'; 'PENN Y'; 'PENN Depth'; 'PENN Radius'; 'PENN DP/mu'};
        case 'SILL'
            nParameters = 7;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'SILL'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.sill{ii}.start;
            model.step(index1:index2) = modelInput.sill{ii}.step;
            model.lower(index1:index2) = modelInput.sill{ii}.lower;
            model.upper(index1:index2) = modelInput.sill{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'SILL Length'; 'SILL Width'; 'SILL Depth'; 'SILL Strike'; 'SILL X'; 'SILL Y'; 'SILL Opening'};
        case 'DIKE'
            nParameters = 8;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'DIKE'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.dike{ii}.start;
            model.step(index1:index2) = modelInput.dike{ii}.step;
            model.lower(index1:index2) = modelInput.dike{ii}.lower;
            model.upper(index1:index2) = modelInput.dike{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'DIKE Length'; 'DIKE Width'; 'DIKE Depth'; 'DIKE Dip'; 'DIKE Strike'; 'DIKE X'; 'DIKE Y'; 'DIKE Opening'};
        case 'FAUL'
            nParameters = 9;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'FAUL'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.fault{ii}.start;
            model.step(index1:index2) = modelInput.fault{ii}.step;
            model.lower(index1:index2) = modelInput.fault{ii}.lower;
            model.upper(index1:index2) = modelInput.fault{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'FAUL Length'; 'FAUL Width'; 'FAUL Depth'; 'FAUL Dip'; 'FAUL Strike'; 'FAUL X'; 'FAUL Y'; 'FAUL StrSlip'; 'FAUL DipSlip'};
            
            % Customised sources start here
        case 'HING'
            nParameters = 11;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'FHIN'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.hing{ii}.start;
            model.step(index1:index2) = modelInput.hing{ii}.step;
            model.lower(index1:index2) = modelInput.hing{ii}.lower;
            model.upper(index1:index2) = modelInput.hing{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
        case 'FHIN'
            nParameters = 11;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'FHIN'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.fhin{ii}.start;
            model.step(index1:index2) = modelInput.fhin{ii}.step;
            model.lower(index1:index2) = modelInput.fhin{ii}.lower;
            model.upper(index1:index2) = modelInput.fhin{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'FAUL Length'; 'FAUL1 Width'; 'FAUL1 Depth'; 'FAUL1 Dip'; 'FAUL Strike'; 'FAUL X'; 'FAUL Y'; 'FAUL StrSlip'; 'FAUL DipSlip'; 'FAUL2 Width'; 'FAUL2 Dip'};
            
        case 'SPLT'
            nParameters = 12;
            index2 = index1 + nParameters - 1;
            [~,ii]=ismember(i,find(strcmp(invpar.model,'SPLT'))); % Find which number geometry that we are on, and take that
            model.m(index1:index2) = modelInput.split{ii}.start;
            model.step(index1:index2) = modelInput.split{ii}.step;
            model.lower(index1:index2) = modelInput.split{ii}.lower;
            model.upper(index1:index2) = modelInput.split{ii}.upper;
            model.gaussPrior(index1:index2) = false(nParameters,1);
            model.parName(index1:index2) = {'SPLT Length'; 'SPLT Base'; 'SPLT Locking'; 'SPLT1 Dip'; 'SPLT Strike'; 'SPLT X'; 'SPLT Y'; 'SPLT1 StrSlip'; 'SPLT1 DipSlip'; 'SPLT2 Dip'; 'SPLT2 StrSlip'; 'SPLT2 DipSlip'};
            
            
            
        otherwise
            error('Invalid model')
    end
    mIx(i+1) = mIx(i) + nParameters;
end

% Add other parameters to invert for (e.g., InSAR constant offset, ramp,
% etc.)
clear index1
index1 = index2+1;
clear index2

nParameters = 0;
insarParName = {};

for i=1:length(insar)
    
    if insar{i}.constOffset == 'y';
        nParameters = nParameters + 1; % Constant offset
        if i == 1
            insarParName = {'InSAR Const.'};
        else
            insarParName = [insarParName, 'InSAR Const.'];
        end
    end
    
    if insar{i}.rampFlag == 'y'
        nParameters = nParameters + 2;  % Linear ramp
        insarParName = [insarParName, 'InSAR X-ramp', 'InSAR Y-ramp';];
    end
end

index2 = index1 + nParameters; % this also adds hyperparameter at end of each model vector

model.m(index1:index2) = zeros(nParameters+1,1);
model.step(index1:index2) = ones(nParameters+1,1)*1e-3;
model.lower(index1:index2) = ones(nParameters+1,1)*-100;
model.upper(index1:index2) = ones(nParameters+1,1)*100;
if ~isempty(insar)
    model.parName(index1:index2-1) = insarParName;
end

% Clear unused values
model.m = model.m(1:index2);
model.step = model.step(1:index2);
model.lower = model.lower(1:index2);
model.upper = model.upper(1:index2);

model.mIx = mIx;
