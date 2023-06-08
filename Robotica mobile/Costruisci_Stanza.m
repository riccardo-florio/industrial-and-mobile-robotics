function [stanza] = Costruisci_Stanza(vertici)
%Costruisci_Stanza costruisce la stanza dato un insieme di vertici
%   Dato l'insieme dei vertici che compongono la stanza e gli ostacoli,
%   restituisce l'insieme di punti che formano i lati degli ostacoli

    stanza=[];
    for i=1:4:length(vertici)
        % lato inferiore
        for j=vertici(i,1):vertici(i+1,1)
            stanza=[stanza; j vertici(i,2)];
        end
        % lato superiore
        for j=vertici(i+2,1):vertici(i+3,1)
            stanza=[stanza; j vertici(i+3,2)];
        end
        % lato sinistro
        for j=vertici(i+1,2):vertici(i+2,2)
            stanza=[stanza; vertici(i,1) j];
        end
        % lato destro
        for j=vertici(i+1,2):vertici(i+2,2)
            stanza=[stanza; vertici(i+3,1) j];
        end
    end
    
end

