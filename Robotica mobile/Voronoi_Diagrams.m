function [path] = Voronoi_Diagrams(start,goal,stanza,ostacoli)
%VORONOI_DIAGRAMS
%   Detailed explanation goes here

    % Plot del diagramma di Voronoi
    figure(); hold on; axis equal; axis([-10 110 -10 110]);
    title('Diagramma di Voronoi'); voronoi(stanza(:,1),stanza(:,2));
    
    % Estrazione dei punti del diagramma
    [x,y] = voronoi(stanza(:,1),stanza(:,2));
    x=x'; y=y';
    
    %% Rimozione dei punti interni agli ostacoli ed esterni alla stanza
    fineLista=false;i=0;
    while ~fineLista
        i=i+1;
        if i>size(x,1)
            fineLista=true;
            break;
        end
        P2=[x(i,1) y(i,1)];
        P1=[x(i,2) y(i,2)];
        for o=5:4:length(ostacoli)
            if (ostacoli(o,1)<P1(1) && P1(1)<ostacoli(o+1,1) && ...
                ostacoli(o,2)<P1(2) && P1(2)<ostacoli(o+2,2)) || ...
                (ostacoli(o,1)<P2(1) && P2(1)<ostacoli(o+1,1) && ...
                ostacoli(o,2)<P2(2) && P2(2)<ostacoli(o+2,2)) || ...
                (P1(1)<ostacoli(1,1) || P1(1)>ostacoli(2,1) || ...
                P1(2)<ostacoli(1,2) || P1(2)>ostacoli(3,2) || ...
                P2(1)<ostacoli(1,1) || P2(1)>ostacoli(2,1) || ...
                P2(2)<ostacoli(1,2) || P2(2)>ostacoli(3,2)) || ...
                (P1(1)==P2(1) && P2(1)==P2(2) && P1(2)==P2(2))

                x(i,:)=[];
                y(i,:)=[];
                i=i-1;
                break
            end
        end
    end

    figure(); plot(x(:,1),y(:,1),"x"); plot(x(:,2),y(:,2),"x");axis('equal')
end

