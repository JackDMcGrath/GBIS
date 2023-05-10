function testGBISPath()

% Test if paths to GBIS are set correctly
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

disp(' ')
which GBISrun
which loadInsarData
which mogi
which PlotInsarWrapped
which generateFinalReport
which runInversion
which local2llh
which fitVariogram
disp(' ')
disp('If you see this message, you are ready to go!')