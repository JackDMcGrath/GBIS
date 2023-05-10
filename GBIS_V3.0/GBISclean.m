function [] = GBISclean(inputDir)

% Delete output directory
%
% Usage:  GBISclean('directoryName')
% =========================================================================
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Markov chain Monte Carlo algorithm incorporating the Metropolis alghoritm
% (e.g., Mosegaard & Tarantola, JGR,(1995).
%
% by Marco Bagnardi and Andrew Hooper (COMET, University of Leeds)
% Email: M.Bagnardi@leeds.ac.uk
% Reference: TBA (Bagnardi and Hooper, in prep.)
% =========================================================================
% Last update: 02/05/2017

%% Check for correct input
if nargin == 0
    error('Directory to remove not specified.')
end

%% Display dialog window to confirm deletion
choice = questdlg(['Do really you want to delete the directory "',inputDir,'"?'], 'Warning!', 'Yes', 'No','Yes');

switch choice
    case 'Yes'
        rmdir(inputDir,'s')
    case 'No'
        return
end