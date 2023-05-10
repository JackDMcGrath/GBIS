function [h,ax,BigAx,patches,pax] = plotmatrix_lower(varargin)
%PLOTMATRIX_LOWER Scatter plot matrix.
%   PLOTMATRIX(X,Y) scatter plots the columns of X against the columns
%   of Y.  If X is P-by-M and Y is P-by-N, PLOTMATRIX will produce a
%   N-by-M matrix of axes. PLOTMATRIX(Y) is the same as PLOTMATRIX(Y,Y)
%   except that the diagonal will be replaced by HIST(Y(:,i)). 
%
%   PLOTMATRIX(...,'LineSpec') uses the given line specification in the
%   string 'LineSpec'; '.' is the default (see PLOT for possibilities).  
%
%   PLOTMATRIX(AX,...) uses AX as the BigAx instead of GCA.
%
%   [H,AX,BigAx,P,PAx] = PLOTMATRIX(...) returns a matrix of handles
%   to the objects created in H, a matrix of handles to the individual
%   subaxes in AX, a handle to big (invisible) axes that frame the
%   subaxes in BigAx, a matrix of handles for the histogram plots in
%   P, and a matrix of handles for invisible axes that control the
%   histogram axes scales in PAx.  BigAx is left as the CurrentAxes so
%   that a subsequent TITLE, XLABEL, or YLABEL will be centered with
%   respect to the matrix of axes.
%
%   Example:
%       x = randn(50,3); y = x*[-1 2 1;2 0 1;1 -2 3;]';
%       plotmatrix(y)

%   Clay M. Thompson 10-3-94
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.19.4.6 $  $Date: 2005/06/21 19:37:51 $

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,3,nargs,'struct'));
if nargs == 3
    labels = varargin{3};
    varargin(3)=[];
    nin=2;
else
nin = nargs;
end

sym = '.'; % Default scatter plot symbol.
dohist = 0;

if ischar(args{nin}),
  sym = args{nin};
  if ~strcmpi(sym,'contour')
      [l,c,m,msg] = colstyle(sym); %#ok
      if ~isempty(msg), error(msg); end %#ok
  end
  nin = nin - 1;
end

if nin==1, % plotmatrix(y)
  rows = size(args{1},2); cols = rows;
  x = args{1}; y = args{1};
  dohist = 1;
elseif nin==2, % plotmatrix(x,y)
  rows = size(args{2},2); cols = size(args{1},2);
  x = args{1}; y = args{2};
else
  error('MATLAB:plotmatrix:InvalidLineSpec',...
        'Invalid marker specification. Type ''help plot''.');
end

% Don't plot anything if either x or y is empty
patches = [];
pax = [];
if isempty(rows) || isempty(cols),
   if nargout>0, h = []; ax = []; BigAx = []; end
   return
end

if ndims(x)>2 || ndims(y)>2,
  error(id('InvalidXYMatrices'),'X and Y must be 2-D.')
end
if size(x,1)~=size(y,1) || size(x,3)~=size(y,3),
  error(id('XYSizeMismatch'),'X and Y must have the same number of rows and pages.');
end

% Create/find BigAx and make it invisible
BigAx = newplot(cax);
fig = ancestor(BigAx,'figure');
hold_state = ishold(BigAx);
set(BigAx,'Visible','off','color','none')

if any(sym=='.'),
  units = get(BigAx,'units');
  set(BigAx,'units','pixels');
  pos = get(BigAx,'Position');
  set(BigAx,'units',units);
  markersize = max(1,min(15,round(15*min(pos(3:4))/max(1,size(x,1))/max(rows,cols))));
else
  markersize = get(0,'defaultlinemarkersize');
end

% Create and plot into axes
ax = zeros(rows,cols);
pos = get(BigAx,'Position');
width = pos(3)/cols;
height = pos(4)/rows;
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];
m = size(y,1);
k = size(y,3);
xlim = nan([rows cols 2]);
ylim = nan([rows cols 2]);
cc=[linspace(1,0,64)',linspace(1,0,64)',ones(64,1)]; % white-blue colormap
cc=[1,1,1;jet(63)];
BigAxHV = get(BigAx,'HandleVisibility');
hh=[];
for i=rows:-1:1
  for j=1:1:i,
    axPos = [pos(1)+(j-1)*width pos(2)+(rows-i)*height ...
             width*(1-space) height*(1-space)];
    findax = findobj(fig,'Type','axes','Position',axPos);
    if isempty(findax),
      ax(i,j) = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',fig);
      set(ax(i,j),'visible','on');
    else
      ax(i,j) = findax(1);
    end
    n_bins=min([round(m/15),100]);
    n_bins2=min([round(sqrt(m/4)),50]);
    if i==rows & dohist==1
      [nn,xx] = hist(reshape(y(:,j,:),[m k]),n_bins);
      patches(i,:) = bar(ax(i,j),xx,nn,1,'b','EdgeColor','none');
      shading flat
      %set(ax(i,j),'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
      %set(histax,'xlim',[xlimmin(1,i)-dx(i) xlimmax(1,i)+dx(i)])
      pax(j) = ax(i,j);  % ax handles for histograms
      axis(ax(i,j),'tight')
    elseif strcmpi (sym,'contour')
      [h,bins] = hist3([reshape(x(:,j,:),[m k]),reshape(y(:,i+1,:),[m k])],[n_bins2,n_bins2]);
      %contour(bins{1},bins{2},h','parent',ax(i,j))';
      %scatter(bins{1},bins{2},20,h','filled','parent',ax(i,j))';
      %imagesc(bins{1},bins{2},h','parent',ax(i,j))';
      h=round(h/max(h(:))*64);
      h=ind2rgb(h',cc);
      image(bins{1},bins{2},h,'parent',ax(i,j))';
      axis xy
      %axis image
      axis tight
      ylim(i,j,:) = get(ax(i,j),'ylim');
    else
      hh(i,j,:) = plot(reshape(x(:,j,:),[m k]), ...
                     reshape(y(:,i+1,:),[m k]),sym,'parent',ax(i,j))';
      set(hh(i,j,:),'markersize',markersize);
      axis tight
      ylim(i,j,:) = get(ax(i,j),'ylim');
    end
    %set(ax(i,j),'xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off')
    set(ax(i,j),'xgrid','off','ygrid','off')
    set(gca,'fontweight','bold','fontsize',8,'YAxisLocation','left')
    xlim(i,j,:) = get(ax(i,j),'xlim');
  end
end

xlimmin = nanmin(xlim(:,:,1),[],1); xlimmax = nanmax(xlim(:,:,2),[],1);
ylimmin = nanmin(ylim(:,:,1),[],2); ylimmax = nanmax(ylim(:,:,2),[],2);

% Try to be smart about axes limits and labels.  Set all the limits of a
% row or column to be the same and inset the tick marks by 10 percent.
inset = .15;

for i=1:rows-dohist,
  set(ax(i,i),'ylim',[ylimmin(i,1) ylimmax(i,1)])
  dy = diff(get(ax(i,i),'ylim'))*inset;
  set(ax(i,ax(i,:)~=0),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
  %set(ax(i,ax(i+dohist,:)~=0),'ylim',[ylimmin(i,1) ylimmax(i,1)])
  l=get(ax(i,1),'YLabel');
  lab = split(string(labels{i+1}));
  l.String=join(lab(2:end));
end

l=get(ax(rows,1),'YLabel');
l.String='Frequency';

dx = zeros(1,cols);
for j=1:cols,
  set(ax(j,j),'xlim',[xlimmin(1,j) xlimmax(1,j)])
  if dohist==1
    dx(j)=0;
  else
    dx(j) = diff(get(ax(j,j),'xlim'))*inset;
  end
  set(ax(ax(:,j)~=0,j),'xlim',[xlimmin(1,j)-dx(j) xlimmax(1,j)+dx(j)])
end

for i=1:rows-1
    set(ax(i,1:i),'xticklabel','')
end
for j=2:cols
    set(ax(j:rows,j),'yticklabel','')
    %set(ax(1:j-1,j),'yticklabel','')
end


set(BigAx,'XTick',get(ax(1,1),'xtick'),'YTick',get(ax(1,1),'ytick'), ...
          'userdata',ax,'tag','PlotMatrixBigAx')

%if dohist, % Put a histogram on the diagonal for plotmatrix(y) case
%  for i=rows:-1:1,
%    histax = axes('Position',get(ax(i,i),'Position'),'HandleVisibility',BigAxHV,'parent',fig);
%    [nn,xx] = hist(reshape(y(:,i,:),[m k]),50);
%    patches(i,:) = bar(histax,xx,nn,'hist');
%    set(histax,'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
%    set(histax,'xlim',[xlimmin(1,i)-dx(i) xlimmax(1,i)+dx(i)])
%    pax(i) = histax;  % ax handles for histograms
%    %ax(i,i)=histax;
%  end
%  patches = patches';
%end

% Make BigAx the CurrentAxes
set(fig,'CurrentAx',BigAx)
if ~hold_state,
   set(fig,'NextPlot','replace')
end

% Also set Title and X/YLabel visibility to on and strings to empty
set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
 'String','','Visible','on')

if nargout~=0,
  h = hh;
end
 
function str=id(str)
str = ['MATLAB:plotmatrix:' str];
