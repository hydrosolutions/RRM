function varargout = setupGUI(varargin)
% SETUPGUI MATLAB code for setupGUI.fig
%      SETUPGUI, by itself, creates a new SETUPGUI or raises the existing
%      singleton*.
%
%      H = SETUPGUI returns the handle to a new SETUPGUI or the handle to
%      the existing singleton*.
%
%      SETUPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPGUI.M with the given input arguments.
%
%      SETUPGUI('Property','Value',...) creates a new SETUPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before setupGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to setupGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setupGUI

% Last Modified by GUIDE v2.5 10-Feb-2015 15:41:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @setupGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @setupGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before setupGUI is made visible.
function setupGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = setupGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


%  1.  GET MODEL NAME
function modelName_Callback(hObject, eventdata, handles)

handles.data.modelName = get(hObject,'String');
handles.modelSetupStatus.String = handles.data.modelName;

guidata(hObject, handles);

function modelName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%  2. GET MODEL DIRECTORY
function modelDirectory_Callback(hObject, eventdata, handles)

% Get main directory
handles.data.modelDir = uigetdir('C:\', 'Select location where to store model directory.');


if ~isfield(handles.data,'modelName')
    handles.data.modelName = datestr(today,'YYYY-MM-DD');
end

handles.cPath = pwd;
mPath = strcat([handles.data.modelDir '/'],handles.data.modelName);
handles.mPath = mPath;

if exist(mPath,'dir')
    cd(mPath)
else
    mkdir(mPath);
    cd(mPath)
end

% CREATE DIRECTORIES
mkdir('./src/')
% data
mkdir('./data/')
mkdir('./data/raw/')
mkdir('./data/raw/ndvi')
mkdir('./data/raw/dem')
mkdir('./data/raw/imomo')
mkdir('./data/processed/')
mkdir('./data/processed/temp/')
mkdir('./data/processed/sub/')
mkdir('./data/processed/db/')

% resources
mkdir('./resources/')
mkdir('./resources/geometry/')
mkdir('./resources/restart/')
mkdir('./resources/samples/')
mkdir('./resources/trainingdata/')
mkdir('./resources/trainingresults/')

% results
mkdir('./results/')
clc
mkdir('./results/nowcast/')
mkdir('./results/forecast/')
mkdir('./results/forecast/Q/')
mkdir('./results/forecast/ET/')
mkdir('./results/forecast/S/')
mkdir('./results/forecast/G/')
% prm
mkdir('./prm/')
msgbox('Model directory successfully generated')
cd(handles.cPath);

set(handles.checkDirectory,'Value',1)


guidata(hObject,handles);


%  3. LOAD SHAPE FILE
function loadShape_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile( {'*.shp','ESRI shape file (*.shp)'},'Pick a Shape file (shp)');
fname=fullfile(pathname,filename);
shpFile = shaperead(fname);

% copyfile([pathname filename(1:end-3) '*'],[handles.data.modelDir '/' handles.data.modelName '/resources/geometry/']);
cd([handles.data.modelDir '/' handles.data.modelName '/resources/geometry/'])
% shpFile = shaperead(filename(1:end-4));
% shpFile = shaperead(shapeFile.shp);

try
    
    mapshow(shpFile);
    xlabel('easting'), ylabel('northing'), grid minor
catch
    errordlg('Error reading shp-files! Please check.')
end

%Write Shapefile
shapewrite(shpFile,'shapeFile.shp');

handles.data.shpFile = shpFile;
bBoxData = [shpFile.BoundingBox];
bBoxR = reshape(bBoxData,[2 2 size(bBoxData,2)/2]);
lon = squeeze(bBoxR(:,1,:)); lon = lon(:);
lat = squeeze(bBoxR(:,2,:)); lat = lat(:);
handles.data.bBox = [min(lon) min(lat) max(lon) max(lat)]; % defined according to the standard bounding box definition


%Round
handles.data.bBox(1)= floor((handles.data.bBox(1)*10))/10;
handles.data.bBox(3)= ceil((handles.data.bBox(3)*10))/10;

if handles.data.bBox(2)<0
    handles.data.bBox(2)= floor((handles.data.bBox(2)*10))/10;
else
    handles.data.bBox(2)= ceil((handles.data.bBox(2)*10))/10;
end

if handles.data.bBox(4)<0
    handles.data.bBox(4)= ceil((handles.data.bBox(4)*10))/10;
else
    handles.data.bBox(4)= floor((handles.data.bBox(4)*10))/10;
end

% Get strings
handles.minLon.String = num2str(handles.data.bBox(1),4);
handles.maxLon.String = num2str(handles.data.bBox(3),4);
handles.minLat.String = num2str(handles.data.bBox(2),4);
handles.maxLat.String = num2str(handles.data.bBox(4),4);

%Update edit box (bounding box)
set(findobj('Tag','maxLon'),'String',handles.maxLon.String)
set(findobj('Tag','minLon'),'String',handles.minLon.String)
set(findobj('Tag','minLat'),'String',handles.minLat.String)
set(findobj('Tag','maxLat'),'String',handles.maxLat.String)


msgbox('Shapefile successfully selected')
cd(handles.cPath);
set(handles.checkShape,'Value',1)

guidata(hObject,handles);


% 4. LOAD CONNECTION MATRIX
function loadConn_Callback(hObject, eventdata, handles)

[filename,pathname] = uigetfile('*.mat','Select a connection matrix');
copyfile([pathname filename(1:end-3) '*'],[handles.data.modelDir '/' handles.data.modelName '/resources/geometry/']);
cd([handles.data.modelDir '/' handles.data.modelName '/resources/geometry/'])
msgbox('Connection matrix successfully selected')
cd(handles.cPath);
set(handles.checkConn,'Value',1)

guidata(hObject,handles);



% 5. UPDATE BOUNDING BOX
function minLat_Callback(hObject, eventdata, handles)

handles.data.bBox(2)=str2num(get(hObject,'String'));

guidata(hObject,handles);
function minLat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxLat_Callback(hObject, eventdata, handles)

handles.data.bBox(4)=str2num(get(hObject,'String'));

guidata(hObject,handles);
function maxLat_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minLon_Callback(hObject, eventdata, handles)

handles.data.bBox(1)=str2num(get(hObject,'String'));

guidata(hObject,handles);
function minLon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxLon_Callback(hObject, eventdata, handles)

handles.data.bBox(3)=str2num(get(hObject,'String'));

guidata(hObject,handles);
function maxLon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% 6. SELECT PRECIPITATION DATA SOURCE
function precip_Callback(hObject, eventdata, handles)
function precip_CreateFcn(hObject, eventdata, handles)

set(hObject,'String',{'FEWS','TRMM'});

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% 7. SELECT TEMPERATURE DATA SOURCE
function temperature_Callback(hObject, eventdata, handles)
function temperature_CreateFcn(hObject, eventdata, handles)

set(hObject,'String',{'GDAS'});

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% 8. SELECT FORECAST DATA SOURCE
function forecast_Callback(hObject, eventdata, handles)
function forecast_CreateFcn(hObject, eventdata, handles)

set(hObject,'String',{'GFS'});

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% 9. SELECT NDVI DATA AREA
function NDVI_Callback(hObject, eventdata, handles)
function NDVI_CreateFcn(hObject, eventdata, handles)

set(hObject,'String',{'N/A','East Africa','North Africa','West Africa','South Africa','Central America','Central Asia','Yemen'});


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% 10. CHECKBOXES
function checkDirectory_Callback(hObject, eventdata, handles)
function checkShape_Callback(hObject, eventdata, handles)
function checkConn_Callback(hObject, eventdata, handles)



% 11. SAVE AND WRITE INFORMATION TO FILE
function saveExit_Callback(hObject, eventdata, handles)


%Write everything in setup
setup.mPath = handles.mPath;
setup.shp = handles.data.shpFile;
setup.bBox = handles.data.bBox;
setup.Path=handles.data.modelDir;

listP=get(handles.precip,'String');
valP=get(handles.precip,'Value');
setup.P=listP{valP};

listT=get(handles.temperature,'String');
valT=get(handles.temperature,'Value');
setup.T=listT{valT};

listF=get(handles.forecast,'String');
valF=get(handles.forecast,'Value');
setup.F=listF{valF};

listN=get(handles.NDVI,'String');
valN=get(handles.NDVI,'Value');
setup.ndvi=listN{valN};

% COPY MAIN FILES
cd(setup.mPath)
cd('../../src/setup')
copyfile('getRaw.m',strcat(setup.mPath,'/src/getRaw_',handles.data.modelName,'.m'))
copyfile('processRaw.m',strcat(setup.mPath,'/src/processRaw_',handles.data.modelName,'.m'))
copyfile('runModel.m',strcat(setup.mPath,'/src/runModel_',handles.data.modelName,'.m'))
copyfile('sendtoDB.m',strcat(setup.mPath,'/src/sendtoDB_',handles.data.modelName,'.m'))

% CREATE START FILE
cd(strcat(setup.mPath,'/src/'))
fid = fopen('master.m', 'w');
fprintf(fid, 'gR=timer(''TimerFcn'',''getRaw_%s'');\n',handles.data.modelName);
fprintf(fid, 'pR=timer(''TimerFcn'',''processRaw_%s'');\n',handles.data.modelName);
fprintf(fid, 'rM=timer(''TimerFcn'',''runModel_%s'');\n\n',handles.data.modelName);
fprintf(fid, 'sD=timer(''TimerFcn'',''sendtoDB_%s'');\n\n',handles.data.modelName);
fprintf(fid, 'set(gR, ''ExecutioMode'', ''fixedRate'') \n');
fprintf(fid, 'set(gR, ''Period'', 10800); \n');
fprintf(fid, 'set(pR, ''ExecutionMode'', ''fixedRate'') \n');
fprintf(fid, 'set(pR, ''Period'', 10800); \n');
fprintf(fid, 'set(rM, ''ExecutionMode'', ''fixedRate'') \n');
fprintf(fid, 'set(rM, ''Period'', 10800); \n\n');
fprintf(fid, 'set(sD, ''Executionode'', ''fixedRate'') \n');
fprintf(fid, 'set(sD, ''Period'', 10800); \n\n');
fprintf(fid, 'startat(gR, ''17:30:00''); \n');
fprintf(fid, 'startat(pR, ''17:40:00''); \n');
fprintf(fid, 'startat(rM, ''17:50:00''); \n');
fprintf(fid, 'startat(sD, ''18:00:00''); \n');
fclose(fid);

% COPY ENKF PARAMETER FILES
cd(setup.mPath)
cd('../../src/enkf/prm')
copyfile('prm-Custom.txt',strcat(setup.mPath,'/prm/prm-Custom.txt'))
copyfile('prm-Global.txt',strcat(setup.mPath,'/prm/prm-Global.txt'))
copyfile('Samples.mat',strcat(setup.mPath,'/resources/samples/samples.mat'))


cd(setup.mPath)
cd('./src')

%save setup to file
save setup.mat setup

close all; 
