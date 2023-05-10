function H=drawmodel(varargin)
%DRAWMODEL     Draws dislocation models
%   H=drawmodel(M) draws the dislocation models specified in
%   the columns of M. Output H is a matrix of handles to the
%   patch and line objects comprising the depicted dislocations.
%
%   H=drawmodel(M,'option',option_val, ...) draws the models
%   using the options specified by pairs of option names and
%   values.  Available options are:
%
%      'color',        [ b | g | {r} | c | m | y | k ]
%      'origin',       [ {'none'} | 2x1 decimal lon,lat ]
%      'outline',      [ {'yes'} | 'no' ]
%      'updipline',    [ {'yes'} | 'no' ]
%      'slipvectors',  [ '2D' | '3D' | {'no'} ]
%      'slipcolor',    [ {b} | g | r | c | m | y | k ]
%      'slipscale',    [ {1} | any real number ]
%      'projection'    [ {3D} | 'no' ]
%      'openingcolor'  [ 'yes' | {'no'} | clim
           
%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 18, 2000  Peter Cervelli        Added more options, standardized code
%   Sep 07, 2000  Peter Cervelli		Original Code
%   Jan 25, 2002  Andy Hooper           Add ability to plot slip vectors
%   Apr 18, 2002  Andy Hooper           Add ability to project on surface
%   Nov 26, 2008  Andy Hooper           Add ability to plot opening as color
%   Dec 05, 2008  Andy Hooper           Add ability to project on surface
%
%-------------------------------------------------------------------------------

%Parse input arguments

    if nargin < 1
        help drawmodel
        return
    end

    if mod((nargin-1),2)
        error('Invalid parameter/value pair arguments.')
    end

    m=varargin{1};
 
%Override defaults with options, if appropriate

    defaults=struct('outline', 'yes', ...
                    'origin', 'none', ...
                    'color', 'r', ...
                    'updipline', 'yes', ...
                    'slipvectors', 'no', ...
                    'slipcolor', 'b',...
                    'slipscale', 1,...
                    'openingcolor', 'no',...
                    'projection','3D');

    opt=optionset(defaults,varargin{2:end});

%Create axes if necessary; otherwise add to current axis

	if isempty(get(gcf,'CurrentAxes'));
	    axes('NextPlot','add');
    	PlotProp='add';
    	view(3)
	else
    	PlotProp=get(gca,'NextPlot');
    	set(gca,'NextPlot','add')
    end

    if ~strcmpi(opt.openingcolor,'no')
        cmap=colormap;
        if isnumeric(opt.openingcolor)
            clim=opt.openingcolor;
        else
            clim=[min(m(10,:)),max(m(10,:))];
        end
        cix=ceil((m(10,:)-clim(1))/(clim(2)-clim(1))*size(cmap,1));
        cix(isnan(cix))=1;
        cix(cix<=0)=1;
        cix(cix>size(cmap,1))=size(cmap,1);
        set(gca,'clim',clim)
    end
        
    
%Loop over models
   
    for i=size(m,2):-1:1

        %Define model corners and slip vector

			vertices=[       0         0   -m(2,i)   -m(2,i)	
                      m(1,i)/2 -m(1,i)/2 -m(1,i)/2  m(1,i)/2
                             0         0         0         0];
                         
            if ~strcmp(opt.slipvectors,'no') 
                slipvec=[-m(9,i);m(8,i);0];  
            else
                slipvec=[0;0;0];
            end
            
        %Create rotation matrices and rotate
               
			sp=sin(m(5,i)*pi/180);
			cp=cos(m(5,i)*pi/180);
			cp(abs(cp)<1e-12)=0;
			sp(abs(sp)<1e-12)=0;
               
			cd=cos(m(4,i)*pi/180);
			sd=sin(m(4,i)*pi/180);
			cd(abs(cd)<1e-12)=0;
			sd(abs(sd)<1e-12)=0;
               
			R2=[cd 0 sd;0 1 0;-sd 0 cd];
			R1=[cp sp 0;-sp cp 0; 0 0 1];
			
			vertices=R1*R2*vertices+repmat([m(6,i) m(7,i) -m(3,i)]',1,4);
            
            slipvec=R1*R2*slipvec;
            slipvectors(:,i)=slipvec;

            centerpoint(:,i)=mean(vertices')';    % store center point of patch
            
     	%Convert to llh if origin is supplied       
        
			if ~strcmp(opt.origin,'none') 
			    vertices(1:2,:)=local2llh(vertices,opt.origin);
                centerpoint(1:2,i)=local2llh(centerpoint(:,i),opt.origin);
                if ~strcmp(opt.slipvectors,'2D')
                    endvec=local2llh(centerpoint(:,i)+slipvec/1000,opt.origin);
                    slipvectors(1:2,i)=endvec-centerpoint(1:2,i);
                end
            end
            
        % 2D and viewing angle
%             if ~strcmp(opt.projection,'3D') 
%                 
%                 az = (opt.projection(1))*pi/180;
%                 el = -opt.projection(2)*pi/180;
%                 normalvec=[-sin(az)*cos(el),cos(az)*cos(el),sin(el)];
%                 if sum(abs(abs(normalvec)-[0,1,0]))<1e-9
%                     yvec=[normalvec(2),0,0];
%                 else
%                     yvec=cross([cos(az),sin(az),0],normalvec)
%                     yvec=yvec/norm(yvec);
%                 end
%                 xvec=cross(normalvec,yvec)
%                 uv = vertices'*[xvec',yvec'];
%                 vertices(1:2,:)=uv';
%             end
    

		%Draw patch and thick line at updip side
		

            if strcmp(opt.outline,'yes')
                facecolor='none';
                edgecolor=opt.color;
            else
                facecolor=opt.color;
                edgecolor='none';
            end
            
            if ~strcmp(opt.openingcolor,'no')
                facecolor=cmap(cix(i),:);
            end
            
            if strcmp(opt.projection,'3D') 
                h(1,i)=patch('Xdata',          vertices(1,:), ...
                             'YData',          vertices(2,:), ...
                             'ZData',          vertices(3,:), ...
                             'HitTest',        'off', ...
                             'FaceColor',      facecolor, ...
                             'EdgeColor',      edgecolor);
            else
                   h(1,i)=patch('Xdata',          vertices(1,:), ...
                             'YData',          vertices(2,:), ...
                             'HitTest',        'off', ...
                             'FaceColor',      facecolor, ...
                             'EdgeColor',      edgecolor);
                   view([0,90])      
            end
                      
                         
            if strcmp(opt.updipline,'yes')
				[Y,I]=sort(vertices(3,:));
                
				if strcmp(opt.projection,'3D') 
				    h(2,i)=plot3(vertices(1,I(3:4)),vertices(2,I(3:4)),vertices(3,I(3:4)), ...
                          'LineWidth',      2, ...
                          'Color',          opt.color);
                else     
                    h(2,i)=plot(vertices(1,I(3:4)),vertices(2,I(3:4)), ...
                          'LineWidth',      2, ...
                          'Color',          opt.color); 
                end
            end

    end
		 
    if ~strcmp(opt.slipvectors,'no') 
         ix=slipvectors(1,:)~=0 & slipvectors(2,:)~=0;
         if strcmp(opt.slipvectors,'2D')
            if ~strcmp(opt.origin,'none')
              quiverll(centerpoint(1,ix),centerpoint(2,ix),slipvectors(1,ix),slipvectors(2,ix),opt.slipscale,0,opt.slipcolor);
            else
              quiverfh(centerpoint(1,ix),centerpoint(2,ix),slipvectors(1,ix),slipvectors(2,ix),opt.slipscale,0,opt.slipcolor);
            end
         else
            quiver3(centerpoint(1,ix),centerpoint(2,ix),centerpoint(3,ix),slipvectors(1,ix)*opt.slipscale,slipvectors(2,ix)*opt.slipscale,slipvectors(3,ix)*opt.slipscale,0,opt.slipcolor)
        end
    end

%Restore plot to previous NextPlot value

    %set(gca,'NextPlot',PlotProp)

%Return handles if called for

	if nargout==1
	    H=h;
	end
