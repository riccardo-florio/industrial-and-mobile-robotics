function [path] = Artificial_Potential_Fields(start,goal,stanza,ostacoli)
%ARTIFICIAL_POTENTIAL_FIELDS
%   Calcola il percorso utilizzando il metodo dei campi potenziali
%   artificiali

    X0 = start;
    G = goal;
    O = stanza;

    %% Calcolo del potenziale totale
    % Definizione delle funzioni Ja, potenziale attrattivo, e 
    % Jr, potenziale repulsivo
    Ja=@(x,y,Gx,Gy)((1/2)*((x-Gx).^2+(y-Gy).^2));
    Jr=@(x,y,Ox,Oy)(1./((x-Ox).^2+(y-Oy).^2));
    
    % Definizione dei gradienti delle funzioni soprastanti
    nablaJaX=@(x,y,Gx,Gy)(x-Gx);
    nablaJaY=@(x,y,Gx,Gy)(y-Gy);
    nablaJrX=@(x,y,Ox,Oy)(2*(Ox-x)./(((x-Ox).^2+(y-Oy).^2)));
    nablaJrY=@(x,y,Ox,Oy)(2*(Oy-y)./(((x-Ox).^2+(y-Oy).^2)));
    
    % Regione di validita'
    dmin=3;
    rho=@(x,y,Ox,Oy)((x-Ox).^2+(y-Oy).^2<=dmin^2);
    
    % Pesi rispettivamente di potenziale attrattivo e repulsivo
    wa=1;
    wo=100;
    
    % Definizione dello spazio di lavoro
    xm=min(X0(1),G(1));xm=min(xm,min(O(:,1)));
    xM=max(X0(1),G(1));xM=max(xM,max(O(:,1)));
    ym=min(X0(2),G(2));ym=min(ym,min(O(:,2)));
    yM=max(X0(2),G(2));yM=max(yM,max(O(:,2)));
    
    deltaXY=1;
    xx=xm:deltaXY:xM;
    yy=ym:deltaXY:yM;
    [XX,YY]=meshgrid(xx,yy);

    % POTENZIALE ATTRATTIVO
    Za=Ja(XX,YY,G(1),G(2));
    nablaJaXX=nablaJaX(XX,YY,G(1),G(2));
    nablaJaYY=nablaJaY(XX,YY,G(1),G(2));
    % Plot
    figure(); grid; surf(XX,YY,Za); title('Potenziale Attrattivo');
    % Normalizzazione di Ja (solo per visualizzarlo meglio)
    nablaJaXXn=nablaJaXX./sqrt(nablaJaXX.^2+nablaJaYY.^2);
    nablaJaYYn=nablaJaYY./sqrt(nablaJaXX.^2+nablaJaYY.^2);
    figure(); grid; quiver(XX,YY,-nablaJaXXn,-nablaJaYYn);
    title('Antigradiente Potenziale Attrattivo'); axis('equal');
    axis([goal(1)-15 goal(1)+15 goal(2)-15 goal(2)+15]);
    
    %POTENZIALE REPULSIVO
    Zr=zeros(size(Za));
    nablaJrXX=zeros(size(nablaJaXX));
    nablaJrYY=zeros(size(nablaJaXX));
    for i=1:size(O,1)
        oi=O(i,:);
        Zr=Zr+Jr(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2));
        nablaJrXX=nablaJrXX+nablaJrX(XX,YY,oi(1),oi(2)).*rho(XX,YY,...
                                                            oi(1),oi(2));
        nablaJrYY=nablaJrYY+nablaJrY(XX,YY,oi(1),oi(2)).*rho(XX,YY,...
                                                            oi(1),oi(2));
    end
    % Si annulla il potenziale all'interno degli ostacoli
    for o=5:4:length(ostacoli)
        for i=1:length(XX)
            for j=1:length(YY)
                if (ostacoli(o,1)<XX(i,j) && XX(i,j)<ostacoli(o+1,1)) ...
                    && (ostacoli(o,2)<YY(i,j) && YY(i,j)<ostacoli(o+2,2))
                    nablaJaXX(i,j)=0;nablaJaYY(i,j)=0;
                    nablaJrXX(i,j)=0;nablaJrYY(i,j)=0;
                end
            end
        end
    end
    % Plot
    figure(); grid; surf(XX,YY,Zr); title('Potenziale Repulsivo');
    % Normalizzazione di Ja (solo per visualizzarlo meglio)
    nablaJrXXn=nablaJrXX./sqrt(nablaJrXX.^2+nablaJrYY.^2);
    nablaJrYYn=nablaJrYY./sqrt(nablaJrXX.^2+nablaJrYY.^2);
    figure(); grid; quiver(XX,YY,-nablaJrXXn,-nablaJrYYn);
    title('Antigradiente Potenziale Repulsivo'); axis('equal');
    axis([0 100 0 100]);
    
    % POTENZIALE TOTALE
    J=wa*Za+wo*Zr;
    nablaJx=wa*nablaJaXX+wo*nablaJrXX;
    nablaJy=wa*nablaJaYY+wo*nablaJrYY;
    
    % Normalizzazione di J (solo per visualizzarlo meglio)
    nablaJxn=nablaJx./sqrt(nablaJx.^2+nablaJy.^2);
    nablaJyn=nablaJy./sqrt(nablaJx.^2+nablaJy.^2);
    
    figure(); grid; surf(XX,YY,J); title('Potenziale Totale');
    figure(); grid; quiver(XX,YY,-nablaJxn,-nablaJyn);
    title('Antigradiente Potenziale Totale'); axis('equal');
    axis([0 100 0 100]);
    
    %% APF
    alpha=0.001;
    TH=0.2;
    nIter=10000;
    path=zeros(2,nIter);
    path(:,1)=X0;
    for k=2:nIter
        Xcur=path(:,k-1);
        nablaA=[
            nablaJaX(Xcur(1),Xcur(2),G(1),G(2));
            nablaJaY(Xcur(1),Xcur(2),G(1),G(2))
            ];
        nablaR=[0;0];
        for i=1:size(O,1)
            oi=O(i,:);
            nablaR=nablaR+[
                nablaJrX(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),...
                                                    Xcur(2),oi(1),oi(2));
                nablaJrY(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),...
                                                    Xcur(2),oi(1),oi(2));
                ];
        end
        nablaJ=wa*nablaA+wo*nablaR;
        Xsucc=Xcur-alpha*nablaJ;
        path(:,k)=Xsucc;
        if norm(Xsucc-G)<=TH
            break;
        end
    end
    if k<nIter
        path(:,k+1:end)=[];
    end
    
    % Plot
    path=path';
    figure(2); hold on; plot(path(:,1),path(:,2),'-'); 
    title('Percorso con Artificial Potential Fields');
end

