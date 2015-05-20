
clear all

tic

%Today's date
today=datenum(date);

% Do the JAVA thing
cd('C:\Users\Sebastian\Dropbox (hydrosolutions)\RRM\src\dbase')
login = 'hydrosolution';
password = 'd8WmFgEGpZ';
url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
namespace= 'http://webservice/';
yesterday=datestr(datenum(date)-1,'yyyy-mm-dd');
javaclasspath('matlabSoapSecurity.jar');
secureSOAP = SecureSOAP(login, password);

%date vector
days=today-6:today;

%Day 1
day1=datestr(datenum(days(1)),'yyyy-mm-dd');
values = secureSOAP.callHTTPSSoap(url, namespace, 'getDateValues',{day1},{'date'}, {'{http://www.w3.org/2001/XMLSchema}string'});

%delete forecast entries
if ischar(values)==0
id1=find([values.variableID]>9);
id2=find([values.variableID]<28);
id3=intersect(id1,id2);
values(id3)=[];
else
values=[];
end
data = values;

%Other days
for m=2:length(days)
day2=datestr(datenum(days(m)),'yyyy-mm-dd');
values = secureSOAP.callHTTPSSoap(url, namespace, 'getDateValues',{day2},{'date'}, {'{http://www.w3.org/2001/XMLSchema}string'});
if ischar(values)==0
%delete forecast entries
id1=find([values.variableID]>9);
id2=find([values.variableID]<28);
id3=intersect(id1,id2);
values(id3)=[];
else
values=[];
end

lv=length(values);
data(end+1:end+lv)=values;
end

% Get site details
for l=1:length(data)
    siteD{l} = secureSOAP.callHTTPSSoap(url, namespace, 'getSiteDetails',{data(l).siteID},{'siteID'}, {'{http://www.w3.org/2001/XMLSchema}int'});
    data(l).siteName=siteD{l}.siteName;
end

% Get user name
load('user.mat')
uID=cell2mat(user(:,1));
uN=user(:,2);
for l=1:length(data)
    name{l}=(uN(uID==data(l).userID));
    data(l).userName=name{l};
end

% Write csv file of data sent
DL(2:length(data)+1,:)=struct2cell(data)';
%empty cells
ec=cellfun(@isempty,DL);
DL(ec)={'N/A'};
DL(1,1)={'Data Value'};
DL(1,2)={'Variable ID'};
DL(1,3)={'Date'};
DL(1,4)={'Latitude'};
DL(1,5)={'Longitude'};
DL(1,6)={'Site ID'};
DL(1,7)={'User ID'};
DL(1,8)={'Site Name'};
DL(1,9)={'User Nickname'};
f=cell2table(DL);
% writetable(f,strcat(num2str(today),'DLoad.csv'))

%load locs
cd('C:\Users\Sebastian\Dropbox (hydrosolutions)\RRM\app\Themi\src')
load('locs.mat')

%Identify locations where no data was sent
for l=1:size(locs,1)
    SID=find([data.siteID]==cell2mat(locs(l,1)));
    countS(l)=length(data(SID));
end
ndID=find(countS==0);

%siteID
ndS=locs(:,1);
ndS=ndS(ndID);
for l=1:length(ndS)
    nd(l).site=cell2mat(ndS(l));
end

%siteName
ndS=locs(:,2);
ndS=ndS(ndID);
for l=1:length(ndS)
    nd(l).siteN=cell2mat(ndS(l));
end

%user
ndS=locs(:,3);
ndS=ndS(ndID);
for l=1:length(ndS)
    nd(l).user=cell2mat(ndS(l));
end

%telephone
ndS=locs(:,4);
ndS=ndS(ndID);
for l=1:length(ndS)
    nd(l).tel=num2str(cell2mat(ndS(l)));
end

% Write csv file of data sent
ND(2:length(nd)+1,:)=(struct2cell(nd'))';
%empty cells
ND(1,1)={'Site ID'};
ND(1,2)={'Site Name'};
ND(1,3)={'User'};
ND(1,4)={'Telephone'};
fND=cell2table(ND);
% writetable(f,strcat(num2str(today),'DLoad.csv'))








%% SEND EMAIL
%Settings
sender = 'imomowb@yahoo.com';
psswd = 'iMoMo567';
recipient={'stoll@hydrosolutions.ch','siegfried@hydrosolutions.ch'};

% YAHOO PREFERENCES 
setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.mail.yahoo.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);
 
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Reduce size for email text
mND=ND;
b=mND(:,2);
c=char(b); 
c(:,20:size(c,2))='';
mND(:,2)=cellstr(c);

mess = evalc('disp(mND)');
subject='No data was sent from the following locations during the last 7 days:';
sendmail(recipient, subject, mess);





toc































