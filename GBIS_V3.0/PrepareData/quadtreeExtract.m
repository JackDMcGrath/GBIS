function [values] = quadtreeExtract(xy,val,xlims,ylims)

% Function to extract values based on quadtree results. Modified after Gonzalez (2015).
%
% Usage: [values]=quadtreeExtract(xy,val,xlims,ylims)
% Input Parameters:
%   xy    : Mx3 array of pixels positions
%   val   : Mx1 array of the value of the pixels
%   xlims : nx2 array of x-coordinate limits for each selected polygon
%   ylims : nx2 array of y-coordinate limits for each selected polygon
%
% Output Parameters:
%   values: n   mean value inside each polygon
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

values    = zeros(1,length(xlims));
  
  % Loop over all quadtree subdivisions
  for i=1:length(xlims)
    in = find( xy(:,2)<=xlims(i,3) & xy(:,2)>=xlims(i,1) & xy(:,3)<=ylims(i,2) & xy(:,3)>=ylims(i,1));
    values(:,i) = mean(val(in));
  end

end