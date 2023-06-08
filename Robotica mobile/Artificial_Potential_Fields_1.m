function [path] = Artificial_Potential_Fields(start, goal, stanza, ...
                                                                ostacoli)
%ARTIFICIAL_POTENTIAL_FIELDS
%   Calcola il percorso utilizzando il metodo dei campi potenziali
%   artificiali

    % Definizione delle funzioni Ja, potenziale attrattivo, e 
    % Jr, potenziale repulsivo
    Ja=@(x,y,Gx,Gy)((1/2)*((x-Gx).^2+(y-Gy).^2));
    Jr=@(x,y,Ox,Oy)(1./((x-Ox).^2+(y-Oy).^2));
    
    % Definizione dei gradienti delle funzioni soprastanti
    nablaJaX=@(x,y,Gx,Gy)(x-Gx);
    nablaJaY=@(x,y,Gx,Gy)(y-Gy);
    nablaJrX=@(x,y,Ox,Oy)(2*(Ox-x)./(((x-Ox).^2+(y-Oy).^2)));
    nablaJrY=@(x,y,Ox,Oy)(2*(Oy-y)./(((x-Ox).^2+(y-Oy).^2)));

    % regione di validita'
    dmin=3;
    rho=@(x,y,Ox,Oy)((x-Ox).^2+(y-Oy).^2<=dmin^2);
    
    wa=1;
    wo=1000;
    
    deltaXY=1;
    
    xm=min(start(1),goal(1));xm=min(xm,min(stanza(:,1)));
    xM=max(start(1),goal(1));xM=max(xM,max(stanza(:,1)));
    ym=min(start(2),goal(2));ym=min(ym,min(stanza(:,2)));
    yM=max(start(2),goal(2));yM=max(yM,max(stanza(:,2)));
    
    xx=xm-2:deltaXY:xM+2;
    yy=ym-2:deltaXY:yM+2;
    [XX,YY]=meshgrid(xx,yy);

    % POTENZIALE ATTRATTIVO
    Za=Ja(XX,YY,goal(1),goal(2));
    nablaJaXX=nablaJaX(XX,YY,goal(1),goal(2));
    nablaJaYY=nablaJaY(XX,YY,goal(1),goal(2));
    
    % Normalizzazione per vedere meglio
    nablaJaXXn=nablaJaXX./sqrt(nablaJaXX.^2+nablaJaYY.^2);
    nablaJaYYn=nablaJaYY./sqrt(nablaJaXX.^2+nablaJaYY.^2);
    
    % Plot potenziale attrattivo
    figure(1);
    grid;surf(XX,YY,Za); title('Potenziale Attrattivo');
    axis([goal(1)-10 goal(1)+10 goal(2)-10 goal(2)+10 0 20]);
    figure(20);title('Antigradiente Potenziale Attrattivo');
    grid;quiver(XX,YY,-nablaJaXXn,-nablaJaYYn);
    axis([goal(1)-10 goal(1)+10 goal(2)-10 goal(2)+10]);
    
    % POTENZIALE REPULSIVO
    Zr=zeros(size(Za));
    nablaJrXX=zeros(length(nablaJaXX),length(nablaJaXX));
    nablaJrYY=zeros(length(nablaJaYY),length(nablaJaYY));

    for i=1:size(stanza,1)
        oi=stanza(i,:);
        Zr=Zr+Jr(XX,YY,oi(1),oi(2)).*rho(XX,YY,oi(1),oi(2));
        nablaJrXX=nablaJrXX+nablaJrX(XX,YY,oi(1),oi(2)).*...
            rho(XX,YY,oi(1),oi(2));
        nablaJrYY=nablaJrYY+nablaJrY(XX,YY,oi(1),oi(2)).*...
            rho(XX,YY,oi(1),oi(2));
    end

    % Normalizzazione per vedere meglio
    nablaJrXXn=nablaJrXX./sqrt(nablaJrXX.^2+nablaJrYY.^2);
    nablaJrYYn=nablaJrYY./sqrt(nablaJrXX.^2+nablaJrYY.^2);
    
    % Plot potenziale repulsivo
    figure(3);
    grid;surf(XX,YY,Zr); title('Parte di Potenziale Repulsivo');
    axis([60 100 60 100 0 20]);
    figure(4);
    grid;quiver(XX,YY,-nablaJrXXn,-nablaJrYYn);
    axis([-5 105 -5 105]);
    title('Antigradiente Potenziale Repulsivo');
    
    % E' necessario annullare il potenziale all'interno degli stanza,
    % altrimenti il rischio è che il robot possa attraversarli 
    for o=5:4:length(ostacoli)
        for i=1:length(XX)
            for j=1:length(YY)
               if(XX(i,j)>ostacoli(o,1) && XX(i,j)<...
                        ostacoli(o+1,1))...
                && (YY(i,j)>ostacoli(o,2) && YY(i,j)<...
                        ostacoli(o+2,2))
                    nablaJrXX(i,j)=0;nablaJrYY(i,j)=0;nablaJaXX(i,j)=0;
                    nablaJaYY(i,j)=0;
                end
            end
        end 
    end
    
    % POTENZIALE TOTALE
    J=wa*Za+wo*Zr;
    nablaJx=wa*nablaJaXX+wo*nablaJrXX;
    nablaJy=wa*nablaJaYY+wo*nablaJrYY;
    
    % Normalizzazione per vedere meglio
    nablaJXn=nablaJx./sqrt(nablaJx.^2+nablaJy.^2);
    nablaJYn=nablaJy./sqrt(nablaJx.^2+nablaJy.^2);
    
    % Plot antigradiente potenziale totale
    figure(5);
    grid;quiver(XX,YY,-nablaJXn,-nablaJYn);
    axis([goal(1)-10 goal(1)+10 goal(2)-10 goal(2)+10]);
    title('Antigradiente Potenziale Totale');
    
    % APF
    alpha=0.001;
    TH=0.2;
    nIter=10000;
    path=zeros(2,nIter);
    path(:,1)=start;
    for k=2:nIter
        Xcur=path(:,k-1);
        nablaA=[
                nablaJaX(Xcur(1),Xcur(2),goal(1),goal(2));
                nablaJaY(Xcur(1),Xcur(2),goal(1),goal(2))
            ];
        nablaR=[0;0];
        for i=1:size(ostacoli,1)
            oi=ostacoli(i,:);
            nablaR=nablaR+[
                nablaJrX(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),Xcur(2),oi(1),oi(2));
                nablaJrY(Xcur(1),Xcur(2),oi(1),oi(2))*rho(Xcur(1),Xcur(2),oi(1),oi(2));
                ];
        end 
        nablaJ=wa*nablaA+wo*nablaR;
        Xsucc=Xcur-alpha*nablaJ;
        
        path(:,k)=Xsucc;
        if norm(Xsucc-goal)<=TH
            break;
        end
    end
    if k<nIter
        path(:,k+1:end)=[];
    end
    
    figure(2); hold on;
    plot(path(1,:),path(2,:),'-'); axis([0 100 0 100]);
    
    

    
end