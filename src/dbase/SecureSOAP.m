classdef SecureSOAP

    properties
        login;
        password;
        passphrase;
    end
    
    methods
        function obj = SecureSOAP (login, password)
            obj.login = login;
            obj.password = password;
            
            % Create passphrase Base64
            encoder = org.apache.commons.codec.binary.Base64();
            passphraseStr = [login, ':', password];
            obj.passphrase = char(encoder.encode(passphraseStr-0))';
        end;
        
        function values = callHTTPSSoap(obj, url, namespace, methodName, paramsValue, paramsName, paramsType)
            import javax.xml.namespace.QName;
            import javax.xml.parsers.DocumentBuilder;
            import javax.xml.parsers.DocumentBuilderFactory;
            import javax.xml.soap.*;

            import org.w3c.dom.Document;
            import org.w3c.dom.Element;
            import org.w3c.dom.NodeList;
            import org.xml.sax.InputSource;

            import java.io.ByteArrayOutputStream;

            import ch.hearc.matlabsoapsecurity.TrustManagerManipulator;

%             disp('Init');
            messageInstance = MessageFactory.newInstance();
            message = messageInstance.createMessage();
            header = message.getSOAPHeader();
            header.detachNode();

            body = message.getSOAPBody();

            soapFactory = SOAPFactory.newInstance();
            bodyName = soapFactory.createName(methodName,'ns2',namespace);
            bodyElement = body.addBodyElement(bodyName);

%             disp('Convert params');
            if length(paramsName) == length(paramsValue)
                for i=1:length(paramsName)
                    for j=1:length(paramsValue{i})
                        symbol = bodyElement.addChildElement(paramsName{i});
                        valueTemp = paramsValue{i};
                        
                        if isa(valueTemp,'char')==0
                            valueTemp = valueTemp(j);
                        end
                        
                        if strcmp(paramsType{i}, '{http://www.w3.org/2001/XMLSchema}double')
                            valueTemp = double(valueTemp);
                        elseif strcmp(paramsType{i}, '{http://www.w3.org/2001/XMLSchema}int')
                            valueTemp = int32(valueTemp);                         
                        elseif strcmp(paramsType{i}, '{http://www.w3.org/2001/XMLSchema}string')
                            valueTemp = char(valueTemp);
                        end

                        symbol.addTextNode(java.lang.String.valueOf(valueTemp));
                    end
                end
            end
%             disp('Convert params end');
            
            message.getMimeHeaders().addHeader('Authorization', ['Basic ',obj.passphrase]); 

            TrustManagerManipulator.allowAllSSL();
            
%             disp('Connection');
            connectionInstance = SOAPConnectionFactory.newInstance();
            connection = connectionInstance.createConnection();
            response = connection.call(message, url);
        
            connection.close();

%             disp('End connection');
            
            byteArrayOS = ByteArrayOutputStream();
            response.writeTo(byteArrayOS);
            reponseXML =  java.lang.String(byteArrayOS.toByteArray());

%             disp(reponseXML);
%             
            values = parseSoapResponse(reponseXML);

%             disp('End');
        end

    end
    
end

