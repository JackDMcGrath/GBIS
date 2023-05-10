%Prepare InSAR data into the mat-file needed by GBIS.
%Eoin Reddin, June 2019
% Need directory with unwrapped, geocoded interferogram. Also need
% parameter files
clc

% Names of master and slaves, path to directory where they are stored
mastername = '20190628';
slavename = '20190704';
filepath = '/scratch/eejdm/ridgecrest/GBIS/';

%% Set up Paths to directory where interferogram and parameter file are stored 
igname = strcat(mastername,'_',slavename,'.diff.unw.geo') ;             %Names of Interferogram and parameter files                               
parname = strcat(mastername,'.rslc.mli.par');
demname = "EQA.dem_par";

% Read in parameter file
fid=fopen(strcat(filepath,parname));                    % Read in MLI Parameter File
mlipar = textscan(fid,'%s%s%s%s');                      % Read in MLI Parameter File
fclose(fid);                                            % Close file 

heading = str2num(cell2mat(mlipar{2}(24)));             % Retrieve Heading from MLI.PAR  
incidence = str2num(cell2mat(mlipar{2}(42)));           % Retrieve incidence angle from MLI.PAR

%Read in Dem file
fid=fopen(strcat(filepath,demname));                    % Read in DEM Parameter File
dempar = textscan(fid,'%s%s%s%s');                      % Read in DEM Parameter File
fclose(fid);                                            % Close file 

%% Assign heading, incidence,width and length from parameter file

len = str2num(cell2mat(dempar{2}(9)));                  % Nrows in interferogram
wid = str2num(cell2mat(dempar{2}(8)));                  % Nlines in interferogram
lat1 = str2num(cell2mat(dempar{2}(10)));                % Latitude at North of Interferogram
lon0 = str2num(cell2mat(dempar{2}(11)));                % Longitude at East of Interferogram
dlat = str2num(cell2mat(dempar{2}(12)));                % Change in Lat per step
dlon = str2num(cell2mat(dempar{2}(13)));                % Change in Lon per step

% Latitude Array
lat0 = lat1-dlat*(len-1);                               % Latitude at south of Interferogram
lat = lat0:dlat:lat1;                                   % Array of Lat values
lat = repmat(lat(1:len),1,wid);                         % Repeat Lat values by n lon values
Lat = reshape(lat,wid*len,1);                           % Reshape to single column
Lat = Lat(1:10:end);                                    % Downsample
Lat(isnan(Lat))=[];                                     % Remove any Nans
        

% Longitude Array
lon1 = lon0+dlon*(wid-1);                               % Longitude at east of Interferogram
lon = lon0:dlon:lon1;                                   % Array of Lat values
lon = repmat(lon(1:wid),len,1);                         % Repeat Lat values by n lon values
Lon = reshape(lon,wid*len,1);                           % Reshape to single column
Lon = Lon(1:10:end);                                    % Downsample
Lon(isnan(Lon))=[];                                     % Remove any Nans

% Heading Array
Heading = ones(size(Lon)).*heading;                     % Heading value times by size of output 

% Incidence Array
Incidence = ones(size(Lon)).*incidence;                 % Incidence value times by size of output (Downsampled Lon) 


%% read in interferogram

fid=fopen(strcat(filepath,igname),'r');                 % Open path to interferogram
[F,g] = fread(fid,'real*4','ieee-be');                  % Read in Interferogram
ifgm=reshape(F,wid,len);                                % Reshape to size
ifgm = ifgm.';                                          % Transpose to shape
ifgm(ifgm ==0) = 0.00001;                               % replace 0 with 1e-5 (required by GBIS)
ifgm(isnan(ifgm)) = 0.00001;                            % replace Nan with 1e-5 (required by GBIS) 
Phase = ifgm(:);                                        % Convert ifgm to column (concatate columns rather than rows)  
Phase = Phase(1:10:end);                                % Downsample to match rest

savename= strcat(mastername,'_',slavename,'.unw.mat');
prompt = 'Would you like to save output for GBIS (*.mat) (y/n): ';
opt = input(prompt,'s');
if opt == 'y'
    disp('Saving output...');
    save(strcat(filename,savename), 'Lon', 'Lat', 'Phase', 'Incidence', 'Heading')
elseif opt == 'n'
    disp('Stack will not be saved');
else
    disp('Invalid input, stack will not be saved');
end

disp('Done!')



 
