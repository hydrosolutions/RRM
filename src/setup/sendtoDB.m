function sendtoDB
try
tic
clear all

warning off

%Today's date
today=datenum(date);

%load setup file
load('setup.mat')

%% Matrix to Vector
cd(strcat(setup.mPath,'/data/processed/db/'))

today=datenum(date);
load(num2str(today));
dummy=ones(1,numel(P)/size(P,3));

% Dates and Site
dMet=repmat(dateMet,numel(P)/size(P,1),1)';
dSim=repmat(dateSim,numel(Q)/size(Q,1),1)';

for k=1:size(P,3)
    for i=1:size(P,2)
        site(:,i,k)=ones(size(P,1),1)*i;
    end
end
idSite=reshape(site,1,numel(P));

% precipitation
Pvec=reshape(P,1,numel(P));
idP(1:size(P,1)*size(P,2))=dummy*13; %mean P =13
idP(size(P,1)*size(P,2)+1:size(P,1)*size(P,2)*2)=dummy*14; %mean P =14
idP(size(P,1)*size(P,2)*2+1:size(P,1)*size(P,2)*3)=dummy*15; %mean P =15

% temperature
Tvec=reshape(T,1,numel(T));
idT(1:size(T,1)*size(T,2))=dummy*10; %mean T =10
idT(size(T,1)*size(T,2)+1:size(T,1)*size(T,2)*2)=dummy*11; %mean T =11
idT(size(T,1)*size(T,2)*2+1:size(T,1)*size(T,2)*3)=dummy*12; %mean T =12

% Discharge
Qvec=reshape(Q,1,numel(Q));
idQ(1:size(Q,1)*size(Q,2))=dummy*25; %mean T =25
idQ(size(Q,1)*size(Q,2)+1:size(Q,1)*size(Q,2)*2)=dummy*26; %mean T =26
idQ(size(Q,1)*size(Q,2)*2+1:size(Q,1)*size(Q,2)*3)=dummy*27; %mean T =27

% soil moisture
Svec=reshape(S,1,numel(S));
idS(1:size(S,1)*size(S,2))=dummy*22; %mean T =22
idS(size(S,1)*size(S,2)+1:size(S,1)*size(S,2)*2)=dummy*23; %mean T =23
idS(size(S,1)*size(S,2)*2+1:size(S,1)*size(S,2)*3)=dummy*24; %mean T =24

% Groundwater
Gvec=reshape(G,1,numel(G));
idG(1:size(G,1)*size(G,2))=dummy*19; %mean T =19
idG(size(G,1)*size(G,2)+1:size(G,1)*size(G,2)*2)=dummy*20; %max T =20
idG(size(G,1)*size(G,2)*2+1:size(G,1)*size(G,2)*3)=dummy*21; %min T =21

% Evapotranspiration
ETvec=reshape(ET,1,numel(ET));
idET(1:size(ET,1)*size(ET,2))=dummy*16; %mean T =16
idET(size(ET,1)*size(ET,2)+1:size(ET,1)*size(ET,2)*2)=dummy*18; %max T =18
idET(size(ET,1)*size(ET,2)*2+1:size(ET,1)*size(ET,2)*3)=dummy*17; %min T =17


%% DO THE WEBSERVICE THING
login = 'hydrosolution';
password = 'd8WmFgEGpZ';
cd(setup.Path)
cd('../src/dbase/')

javaclasspath('matlabSoapSecurity.jar');
secureSOAP = SecureSOAP(login, password);
url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
namespace= 'http://webservice/';



%% SEND DATA
% precipitation
for l=1:length(Pvec)
    values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValue',{Pvec(l),dMet(l),idSite(l),idP(l)},{'dataValue','dateTimeUTC','siteID','variableID'},{'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});
end
 
% Temperature
for l=1:length(Tvec)
    values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValue',{Tvec(l),dMet(l),idSite(l),idT(l)},{'dataValue','dateTimeUTC','siteID','variableID'},{'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});
end

% Discharge
for l=1:length(Qvec)
    values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValue',{Qvec(l),dSim(l),idSite(l),idQ(l)},{'dataValue','dateTimeUTC','siteID','variableID'},{'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});
end

% Soil moisture
for l=1:length(Svec)
    values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValue',{Svec(l),dSim(l),idSite(l),idS(l)},{'dataValue','dateTimeUTC','siteID','variableID'},{'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});
end

% Groundwater
for l=1:length(Gvec)
    values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValue',{Gvec(l),dSim(l),idSite(l),idG(l)},{'dataValue','dateTimeUTC','siteID','variableID'},{'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});
end
 
% Evapotranspiration
for l=1:length(ETvec)
    values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValue',{ETvec(l),dSim(l),idSite(l),idET(l)},{'dataValue','dateTimeUTC','siteID','variableID'},{'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});
end


%% Create Geotiff and send to database
% get data
cd(strcat(setup.mPath,'/data/raw/'))

load GDAS.mat;
load FEWS.mat;


% create temporary directory and change dir
tempDir = 'temp'; % temporary directory for data download
removeTemp = 1;   % remove temporary directory: yes = 1, no = 0;
mkdir(tempDir);
cd(tempDir);

% create geotiff
geotiffwrite(strcat('P',num2str(timeP(end))),P(:,:,end),R_P);
geotiffwrite(strcat('meanT',num2str(timeT(end))),T(:,:,end),R_T);
geotiffwrite(strcat('maxT',num2str(timeT(end))),T(:,:,end),R_T);
geotiffwrite(strcat('minT',num2str(timeT(end))),T(:,:,end),R_T);

% do the webservice thing
login = 'hydrosolution';
password = 'd8WmFgEGpZ';
cd(setup.Path)
cd('../src/dbase/')
javaclasspath('matlabSoapSecurity.jar');
secureSOAP = SecureSOAP(login, password);
url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
namespace= 'http://webservice/';


% Define files to upload
filepath = strcat(setup.mPath,'/data/raw/temp/');
dvec=datevec(timeT(end));
dvec(4)=12;
datetime=datestr(dvec,'yyyy-mm-dd HH:MM:SS');


% mean T
filename = strcat('meanT',num2str(timeT(end)),'.tif');
variableId = 63;
fid = fopen(strcat(filepath, '/', filename));
fseek(fid, 0,'eof');
filelength = ftell(fid);
fseek(fid,0,'bof');
dataTemp = fread(fid,filelength, '*uchar');
data = transpose(dataTemp);
fclose(fid);

% Encoding data (binary to base64 string)
encoder = org.apache.commons.codec.binary.Base64;
dataBase64Temp = encoder.encodeBase64(data);
dataBase64 = transpose(dataBase64Temp);
dataBase64String = char(dataBase64);
errorCode = secureSOAP.callHTTPSSoap(url, namespace, 'sendSatelliteData',{dataBase64String, filename, datetime, variableId},{'data','filename', 'datetime', 'variableId'},{'{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}int'});

% max T
filename = strcat('maxT',num2str(timeT(end)),'.tif');
variableId = 65;
fid = fopen(strcat(filepath, '/', filename));
fseek(fid, 0,'eof');
filelength = ftell(fid);
fseek(fid,0,'bof');
dataTemp = fread(fid,filelength, '*uchar');
data = transpose(dataTemp);
fclose(fid);

% Encoding data (binary to base64 string)
encoder = org.apache.commons.codec.binary.Base64;
dataBase64Temp = encoder.encodeBase64(data);
dataBase64 = transpose(dataBase64Temp);
dataBase64String = char(dataBase64);
errorCode = secureSOAP.callHTTPSSoap(url, namespace, 'sendSatelliteData',{dataBase64String, filename, datetime, variableId},{'data','filename', 'datetime', 'variableId'},{'{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}int'});

% min T
filename = strcat('minT',num2str(timeT(end)),'.tif');
variableId = 64;
fid = fopen(strcat(filepath, '/', filename));
fseek(fid, 0,'eof');
filelength = ftell(fid);
fseek(fid,0,'bof');
dataTemp = fread(fid,filelength, '*uchar');
data = transpose(dataTemp);
fclose(fid);

% Encoding data (binary to base64 string)
encoder = org.apache.commons.codec.binary.Base64;
dataBase64Temp = encoder.encodeBase64(data);
dataBase64 = transpose(dataBase64Temp);
dataBase64String = char(dataBase64);
errorCode = secureSOAP.callHTTPSSoap(url, namespace, 'sendSatelliteData',{dataBase64String, filename, datetime, variableId},{'data','filename', 'datetime', 'variableId'},{'{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}int'});


% PRECIPITATION
% Define files to upload
dvec=datevec(timeP(end));
dvec(4)=12;
datetime=datestr(dvec,'yyyy-mm-dd HH:MM:SS');

filename = strcat('P',num2str(timeP(end)),'.tif');
variableId = 62;
fid = fopen(strcat(filepath, '/', filename));
fseek(fid, 0,'eof');
filelength = ftell(fid);
fseek(fid,0,'bof');
dataTemp = fread(fid,filelength, '*uchar');
data = transpose(dataTemp);
fclose(fid);

% Encoding data (binary to base64 string)
encoder = org.apache.commons.codec.binary.Base64;
dataBase64Temp = encoder.encodeBase64(data);
dataBase64 = transpose(dataBase64Temp);
dataBase64String = char(dataBase64);
errorCode = secureSOAP.callHTTPSSoap(url, namespace, 'sendSatelliteData',{dataBase64String, filename, datetime, variableId},{'data','filename', 'datetime', 'variableId'},{'{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}int'});


% delete temporary directory
cd(strcat(setup.mPath,'/data/raw/'))
if removeTemp == 1
     rmdir(tempDir, 's')
end
toc
clear all
catch
disp('Error occured, next try in 3 hours')  
clear all
end








    