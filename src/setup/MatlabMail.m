try

clear all
warning off

%Today's date
today=datenum(date);

%load setup file
load('setup.mat')    
    
    
%% Data Status
cd(strcat(setup.mPath,'/data/raw/'))

today=datenum(date);
%FEWS
load('FEWS.mat','timeP');
diffFEWS=today-timeP(end);

%GDAS
load('GDAS.mat','timeT');
diffGDAS=today-timeT(end);

%GFS
load('GFS.mat');
diffGFS=today-F{end}.timeF(1);


message=['iMoMo water balance tool is running' 10 ...
        strcat('Date of the last FEWS data:',datestr(today-diffFEWS)) 10 ...
        strcat('Date of the last GDAS data:',datestr(today-diffGDAS)) 10 ...
        strcat('Date of the last GFS data:',datestr(today-diffGFS))];

subject='Status of iMoMo water balance tool';

% MAIL SETTINGS
sender = 'imomowb@yahoo.com';
psswd = 'iMoMo567';
recipient={'stoll@hydrosolutions.ch'};%;'siegfried@hydrosolutions.ch'};

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


% SEND MAIL
sendmail(recipient, subject, message);




%% Daily Data 
cd(strcat(setup.mPath,'/data/raw/imomo/'))
load(strcat(num2str(today),'DLoad'),'data');

%Settings
sender = 'imomowb@yahoo.com';
psswd = 'iMoMo567';
recipient={'stoll@hydrosolutions.ch','luethi@photrack.ch','rnaudascher@yahoo.de','bambabakary@yahoo.com','siegfried@hydrosolutions.ch'};

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

% Message
if isempty(data)==0
    message='Please find attached the lastet uploads to the iMoMo database!';
    subject='Newest entries to iMoMo database';
    attachment={strcat(num2str(today),'DLoad.csv')};
    sendmail(recipient, subject, message, attachment);
else
    message='Today no data was uploaded to the iMoMo database!';
    subject='Newest entries to iMoMo database';
    sendmail(recipient, subject, message);
end
catch
end


