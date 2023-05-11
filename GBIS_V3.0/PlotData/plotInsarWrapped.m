function plotInsarWrapped(xy, los, wavelength, cmap, name)

% Function to plot wrapped InSAR data
%
% Usage: plotInsarWrapped(xy, los, cmap, name)
% Input Parameters:
%       xy: local coordinates of data points
%       los: line-of-sight displacement measured at data points
%       wavelength: data wavelength for wrapping
%       cmap: colormaps for plotting
%       name: name of dataset for figure title
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

%% Patch scattered data for faster plotting
    edge = round(min(abs(diff(xy(:,3)))))+2; % Size of patch set to minumum distance between points
    if edge < 50
        edge = 50;
    end
    xs = [xy(:,2)'; xy(:,2)'+edge; xy(:,2)'+edge; xy(:,2)']; % Coordinates of four vertex of patch
    ys = [xy(:,3)'; xy(:,3)'; xy(:,3)'+edge; xy(:,3)'+edge];
     
   
    % Display wrapped interferogram at 5.6 cm wavelength
    colormap(cmap.seismo); 
    h1 = patch(xs, ys, 'r');
    set(h1, 'facevertexcdata', mod(los,wavelength/2), 'facecolor', 'flat', 'edgecolor', 'none')
    axis equal; axis tight;
    ax = gca;
    grid on
    ax.Layer = 'top';
    ax.Box = 'on';
    ax.LineWidth = 1.0;
    ax.GridLineStyle = '--';
    cbar = colorbar; ylabel(cbar,'LOS displacement m','FontSize', 14); 
    xlabel('X distance from origin (m)','FontSize', 14)
    ylabel('Y distance from origin (m)','FontSize', 14)
    t = title(['Wrapped ', name],'FontSize', 18);
    set(t,'Position',get(t,'Position')+[0 1000 0]);
    drawnow


    

