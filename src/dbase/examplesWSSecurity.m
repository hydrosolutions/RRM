login = 'hydrosolution';
password = 'd8WmFgEGpZ';

cd('~/Dropbox/iMoMoPHASE1/model/iMoMoToolbox/dbase')
javaclasspath('matlabSoapSecurity.jar');

secureSOAP = SecureSOAP(login, password);

url = 'https://157.26.64.27/HydrosolutionWS_V2/HydrosolutionWS';
namespace= 'http://webservice/';
    
% values = secureSOAP.callHTTPSSoap(url, namespace, 'getDailyValues',{},{}, {});

% values = secureSOAP.callHTTPSSoap(url, namespace, 'specificValuesRequest',{1,'2012-03-18 23:30:59','2013-04-20 23:30:59',-4,-2,20,40},{'variableID','start_dateTimeUTC','end_dateTimeUTC','min_latitude','max_latitude','min_longitude','max_longitude'}, {'{http://www.w3.org/2001/XMLSchema}int', '{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}string', '{http://www.w3.org/2001/XMLSchema}double', '{http://www.w3.org/2001/XMLSchema}double', '{http://www.w3.org/2001/XMLSchema}double', '{http://www.w3.org/2001/XMLSchema}double'});

% values = secureSOAP.callHTTPSSoap(url, namespace, 'getVariableDetails',{1},{'variableID'}, {'{http://www.w3.org/2001/XMLSchema}int'});

% values = secureSOAP.callHTTPSSoap(url, namespace, 'getSiteDetails',{1},{'siteID'}, {'{http://www.w3.org/2001/XMLSchema}int'});

% values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValue',{14.3,'2012-04-04 00:30:00',1,1},{'dataValue','dateTimeUTC','siteID','variableID'}, {'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});

dataValue=[13 22 31];
dateTime={'2013-04-18 10:00:00' '2013-04-18 11:00:00' '2013-04-18 12:00:00'};
siteID=[1 1 1];
variableID=[1 1 1];
values = secureSOAP.callHTTPSSoap(url, namespace, 'sendValueVector',{Pvec,dMet,idSite,idP},{'dataValue','dateTimeUTC','siteID','variableID'},{'{http://www.w3.org/2001/XMLSchema}double','{http://www.w3.org/2001/XMLSchema}string','{http://www.w3.org/2001/XMLSchema}int','{http://www.w3.org/2001/XMLSchema}int'});



