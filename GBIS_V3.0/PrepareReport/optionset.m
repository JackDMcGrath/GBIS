function opt=optionset(varargin)
%OPTIONSET    opt=optionset(defaults,options)
%
%Sets options structure to defaults and overrides if appropriate

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date           Programmer            Description of Change
%   ====           ==========            =====================
%
%   Aug 14, 2001   Peter Cervelli                Original Code
%
%-------------------------------------------------------------

defaults=varargin{1};

if nargin > 1
    options=struct(varargin{2:end});
    optfields=fieldnames(options);
    for i=1:length(optfields)
        if isfield(defaults,lower(optfields{i}))
            defaults=setfield(defaults,lower(optfields{i}),getfield(options,optfields{i}));
        else
            warning(['Unknown option: ',optfields{i}]);
        end
    end
end

opt=defaults;
