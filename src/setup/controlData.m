function controlData
clear all

%Today's date
today=datenum(date);
%load setup file
load('setup.mat') 
cd(setup.Path)
cd('../src/dbase/')


%% To the JAVA thing
login = 'hydrosolution';
password = 'd8WmFgEGpZ';
url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
namespace= 'http://webservice/';

javaclasspath('matlabSoapSecurity.jar');
secureSOAP = SecureSOAP(login, password);


%% Get daily values
values = secureSOAP.callHTTPSSoap(url, namespace, 'getDailyValues',{},{}, {});

% Filter values
id1=find([values.siteID]==25); %Olungura Siwandeti - Lembros 0755816128
count1=length(values(id1));
id2=find([values.siteID]==22); %Olangit Olmringa\Josho - Osambi 0763469951
count2=length(values(id2));
id3=find([values.siteID]==23); %Saitabau\Saumu - Osambi 0763469951
count3=length(values(id3));
id4=find([values.siteID]==19); %Sombetini - Adamson 0752932775
count4=length(values(id4));
id5=find([values.siteID]==20); %Olasiti Kati -  Adamson 0752932776
count5=length(values(id5));
id6=find([values.siteID]==36); %Moivo Mollel - Selemani 0756541574
count6=length(values(id6));
id7=find([values.siteID]==16); %Kichangani - Jumbe Lomitu 0753869298
count7=length(values(id7));
id8=find([values.siteID]==15); %Mungushi - Jackson 0765947605
count8=length(values(id8));
id9=find([values.siteID]==26); %Themi ya Simba - Musa 0768940056
count9=length(values(id9));
id10=find([values.siteID]==21); %Mfereji wa Wazee - Stefano 0752291673
count10=length(values(id10));
id11=find([values.siteID]==42); 
count11=length(values(id11));
id12=find([values.siteID]==45); 
count12=length(values(id12));
id13=find([values.siteID]==27); 
count13=length(values(id13));
id14=find([values.siteID]==41); 
count14=length(values(id14));
id15=find([values.siteID]==43); 
count15=length(values(id15));
id16=find([values.siteID]==40); 
count16=length(values(id16));

%% Load old data
cd(strcat(setup.mPath,'/data/processed/db/'))
delimiter = ',';
startRow = 2;
filename='report.csv';
formatSpec = '%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

for l=1:length(dataArray)
    data(:,l)=dataArray{l};
end

%Data
data{1,end+1}=datestr(today);
data{2,end}=count1;
data{3,end}=count2;
data{4,end}=count3;
data{5,end}=count4;
data{6,end}=count5;
data{7,end}=count6;
data{8,end}=count7;
data{9,end}=count8;
data{10,end}=count9;
data{11,end}=count10;
data{12,end}=count11;
data{13,end}=count12;
data{14,end}=count13;
data{15,end}=count14;
data{16,end}=count15;
data{17,end}=count16;


%% Write csv
%Details
data{1,1}='Site ID';data{1,2}='Person';data{1,3}='Telephone';data{1,4}='Location';
data{2,1}='25';data{2,2}='Lembros';data{2,3}='0755816128';data{2,4}='Olungura Siwandeti';
data{3,1}='22';data{3,2}='Osambi';data{3,3}='0763469951';data{3,4}='Olangit Olmringa\Josho';
data{4,1}='23';data{4,2}='Osambi';data{4,3}='0763469951';data{4,4}='Saitabau\Saumu';
data{5,1}='19';data{5,2}='Adamson';data{5,3}='0752932776';data{5,4}='Sombetini';
data{6,1}='20';data{6,2}='Selemani';data{6,3}='0756541574';data{6,4}='Olasiti Kati';
data{7,1}='36';data{7,2}='Jumbe Lomitu';data{7,3}='0753869298';data{7,4}='Moivo Mollel';
data{8,1}='16';data{8,2}='Thomas';data{8,3}='255759480365';data{8,4}='Kichangani Furrow';
data{9,1}='15';data{9,2}='Julius';data{9,3}='255768582513';data{9,4}='Mungushi Furrow';
data{10,1}='26';data{10,2}='Mussa';data{10,3}='255765826995';data{10,4}='Themi ya Simba';
data{11,1}='21';data{11,2}='Stefano';data{11,3}='0752291673';data{11,4}='Mfereji wa Wazee';
data{12,1}='42';data{12,2}='Mussa';data{12,3}='255765826995';data{12,4}='Downstr. T. ya Simba 1';
data{13,1}='45';data{13,2}='Mussa';data{13,3}='255765826995';data{13,4}='Fili estate Furrow';
data{14,1}='27';data{14,2}='Mussa';data{14,3}='255765826995';data{14,4}='Downstr. T. ya Simba 2';
data{15,1}='41';data{15,2}='Thomas';data{15,3}='255759480365';data{15,4}='Kichangani River';
data{16,1}='43';data{16,2}='Thomas';data{16,3}='255759480365';data{16,4}='Kigongoni Furrow';
data{17,1}='40';data{17,2}='Thomas';data{17,3}='255759480365';data{17,4}='Mungushi River';


cd(strcat(setup.mPath,'/data/processed/db/'))
d=cell2table(data);
writetable(d,'report.csv')

%Settings
sender = 'imomowb@yahoo.com';
psswd = 'iMoMo567';
recipient={'stoll@hydrosolutions.ch','rnaudascher@yahoo.de','stoll@hydrosolutions.ch','siegfried@hydrosolutions.ch',...
    'tjkitomari@yahoo.com','bambabakary@yahoo.com','hoseasanga@gmail.com'};

% YAHOO PREFERENCES 
setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.mail.yahoo.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);
 
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
% props.put('mail.smtp.starttls.enable','true');

message1=['iMoMo water balance tool is running' 10 ...
        strcat('Date of the last FEWS data:',datestr(date)) 10 ...
        strcat('Date of the last GDAS data:',datestr(date)) 10 ...
        strcat('Date of the last GFS data:',datestr(date))];
     

% Message
% message=['Yesterday the following data was submitted:' 10 ...
%         'Site ID' ' ' data{2,1} ': ' data{2,2} ' ' 'at' ' ' data{2,4} ' ' 'sent' ' ' num2str(data{2,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{3,1} ': ' data{3,2} ' ' 'at' ' ' data{3,4} ' ' 'sent' ' ' num2str(data{3,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{4,1} ': ' data{4,2} ' ' 'at' ' ' data{4,4} ' ' 'sent' ' ' num2str(data{4,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{5,1} ': ' data{5,2} ' ' 'at' ' ' data{5,4} ' ' 'sent' ' ' num2str(data{5,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{6,1} ': ' data{6,2} ' ' 'at' ' ' data{6,4} ' ' 'sent' ' ' num2str(data{6,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{7,1} ': ' data{7,2} ' ' 'at' ' ' data{7,4} ' ' 'sent' ' ' num2str(data{7,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{8,1} ': ' data{8,2} ' ' 'at' ' ' data{8,4} ' ' 'sent' ' ' num2str(data{8,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{9,1} ': ' data{9,2} ' ' 'at' ' ' data{9,4} ' ' 'sent' ' ' num2str(data{9,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{10,1} ': ' data{10,2} ' ' 'at' ' ' data{10,4} ' ' 'sent' ' ' num2str(data{10,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{11,1} ': ' data{11,2} ' ' 'at' ' ' data{11,4} ' ' 'sent' ' ' num2str(data{11,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{12,1} ': ' data{12,2} ' ' 'at' ' ' data{12,4} ' ' 'sent' ' ' num2str(data{12,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{13,1} ': ' data{13,2} ' ' 'at' ' ' data{13,4} ' ' 'sent' ' ' num2str(data{13,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{14,1} ': ' data{14,2} ' ' 'at' ' ' data{14,4} ' ' 'sent' ' ' num2str(data{14,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{15,1} ': ' data{15,2} ' ' 'at' ' ' data{15,4} ' ' 'sent' ' ' num2str(data{15,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{16,1} ': ' data{16,2} ' ' 'at' ' ' data{16,4} ' ' 'sent' ' ' num2str(data{16,end}) ' ' 'values' 10 ...
%         'Site ID' ' ' data{17,1} ': ' data{17,2} ' ' 'at' ' ' data{17,4} ' ' 'sent' ' ' num2str(data{17,end}) ' ' 'values'];
    
message=['Yesterday the following data was submitted:' 10 ...
        'Site ID' ' ' data{2,1} ': ' data{2,2} ' ' 'at' ' ' data{2,4} ' ' 'sent' ' ' num2str(data{2,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{3,1} ': ' data{3,2} ' ' 'at' ' ' data{3,4} ' ' 'sent' ' ' num2str(data{3,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{4,1} ': ' data{4,2} ' ' 'at' ' ' data{4,4} ' ' 'sent' ' ' num2str(data{4,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{5,1} ': ' data{5,2} ' ' 'at' ' ' data{5,4} ' ' 'sent' ' ' num2str(data{5,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{6,1} ': ' data{6,2} ' ' 'at' ' ' data{6,4} ' ' 'sent' ' ' num2str(data{6,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{7,1} ': ' data{7,2} ' ' 'at' ' ' data{7,4} ' ' 'sent' ' ' num2str(data{7,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{8,1} ': ' data{8,2} ' ' 'at' ' ' data{8,4} ' ' 'sent' ' ' num2str(data{8,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{9,1} ': ' data{9,2} ' ' 'at' ' ' data{9,4} ' ' 'sent' ' ' num2str(data{9,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{10,1} ': ' data{10,2} ' ' 'at' ' ' data{10,4} ' ' 'sent' ' ' num2str(data{10,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{11,1} ': ' data{11,2} ' ' 'at' ' ' data{11,4} ' ' 'sent' ' ' num2str(data{11,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{12,1} ': ' data{12,2} ' ' 'at' ' ' data{12,4} ' ' 'sent' ' ' num2str(data{12,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{13,1} ': ' data{13,2} ' ' 'at' ' ' data{13,4} ' ' 'sent' ' ' num2str(data{13,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{14,1} ': ' data{14,2} ' ' 'at' ' ' data{14,4} ' ' 'sent' ' ' num2str(data{14,end}) ' ' 'values' 10 ...
        'Site ID' ' ' data{16,1} ': ' data{16,2} ' ' 'at' ' ' data{16,4} ' ' 'sent' ' ' num2str(data{16,end}) ' ' 'values'];    
 
    
subject='Daily report on data sending activity';
attachment={strcat('report.csv')};
sendmail(recipient, subject, message, attachment);


end







