function [cb2,ax1,ax2] = newcolorbar(varargin) 
% The newcolorbar function creates a new set of invisible axes matching
% the size and extents of the current axes.  This allows a second colormap
% to be used in such a way that it is perceived to be two separate colormaps 
% and two colorbars in the same axes. 
% 
%% Syntax 
% 
%  colorbar
%  colorbar(Location)
%  colorbar(...,PropertyName,PropertyValue)
%  [cb2,ax2,ax1] = colorbar(...)
% 
%% Description 
% 
% colorbar creates a new set of axes and a new colorbar in the default 
% (right) location. 
% 
% colorbar(Location) specifies location of the new colorbar as 
%       'North'              inside plot box near top
%       'South'              inside bottom
%       'East'               inside right
%       'West'               inside left
%       'NorthOutside'       outside plot box near top
%       'SouthOutside'       outside bottom
%       'EastOutside'        outside right (default)  
%       'WestOutside'        outside left
%
% colorbar(...,PropertyName,PropertyValue) specifies additional
% name/value pairs for colorbar. 
% 
% [cb2,ax2,ax1] = colorbar(...) returns handles of the new colorbar
% cb2, the new axes ax2, and the previous current axes ax1. 
% 
%% Example 1 
% Let's plot some parula-colored scattered data atop a grayscale gridded
% dataset.  Start by plotting the gridded data.  We'll use the inbuilt
% peaks dataset for this example: 
% 
% pcolor(peaks(500))
% shading interp
% colormap(gca,gray(256))
% colorbar('southoutside')
% 
% The newcolorbar function differs from Matlab's colorbar function in
% that we have to create a newcolorbar before plotting the new
% color-scaled dataset.  We can create a new colorbar quite simply:
% 
% newcolorbar
% 
% Then plot some random color-scaled scattered data: 
% Plot a scattered dataset with parula colormap: 
% 
% scatter(500*rand(30,1),500*rand(30,1),60,100*rand(30,1),'filled') 
% 
%% Example 2 
% Now let's repeat Example 1, but add a little formatting.  I'm using 
% Stephen Cobeldick's brewermap function to create the lovely colormaps: 
% 
% figure
% pcolor(peaks(500))
% shading interp
% colormap(gca,brewermap(256,'greens'))
% cb1 = colorbar('southoutside'); 
% xlabel(cb1,'colorbar for peaks data') 
% 
% Add a new colorbar toward the bottom inside of the current axes and 
% specify blue text: 
% 
% cb2 = newcolorbar('south','color','blue'); 
% 
% Plot scattered data and label the new colorbar: 
% 
% scatter(500*rand(30,1),500*rand(30,1),100,8*randn(30,1),'filled') 
% colormap(gca,brewermap(256,'*RdBu'))
% caxis([-10 10]) % sets scattered colorbar axis 
% xlabel(cb2,'colorbar for scattered data','color','blue')
% 
%% Author Info
% The newcolorbar function was written by  Chad A. Greene of the University
% of Texas at Austin's Institute for Geophysics (UTIG), August 2015. 
% http://www.chadagreene.com. 
% 
% See also colorbar and colormap. 


%% Make sure user has R2014b or later: 

assert(verLessThan('matlab','8.4.0')==0,'Sorry the newcolorbar function requires Matlab R2014b or later.') 

%% Begin work: 

% Axis 1 is the original current axis: 
ax1 = gca; 

% Get position of axis 1: 
ax1p = get(ax1,'pos'); 

% Create new axes: 
ax2 = axes;

% Co-locate ax2 atop ax1: 
set(ax2,'pos',ax1p)

% Make ax2 invisible:
axis off 

% Link ax1 and ax2 so zooming will work properly: 
linkaxes([ax1,ax2],'xy') 

% Create a new colorbar
cb2 = colorbar(varargin{:}); 

% If creation of the new colorbar resized ax2, set ax1 to the same size. 
% But first we need to make sure everything is drawn properly: 
drawnow 
set(ax1,'pos',get(ax2,'pos'))

% Make new axes current: 
axes(ax2)

% Ensure we're ready to plot new data:
hold on

%% Clean up: 

if nargout==0
    clear cb2
end
