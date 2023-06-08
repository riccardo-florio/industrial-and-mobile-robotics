function [path] = Discrete_Potential_Fields(start, goal, ostacoli)
%Discrete_Potential_Fields
%   Calcola il percorso utilizzando il metodo dei campi potenziali discreti

    n = ostacoli(2,1)-ostacoli(1,1); % Numero di righe della matrice
    m = ostacoli(3,2)-ostacoli(2,2); % Numero di colonne della matrice
    xG = goal(1); % Coordinata x del punto di goal
    yG = goal(2); % Coordinata y del punto di goal

    %% APPROSSIMAZIONE DISCRETA DEL POTENZIALE
    % Matrice con tutti i valori pari a -1
    griglia = zeros(n,m);
    for i=1:n
        for j=1:m
            griglia(i,j)=-1;
        end
    end

    % 0 nel goal
    griglia(xG,yG)=0;

    for k=1:m*n
        for i=1:n
            for j=1:m
                if griglia(i,j)==k-1
                    % nw
                    if i-1>0 && j-1>0 && griglia(i-1,j-1)==-1
                        griglia(i-1,j-1)=k;
                    end
                    % n
                    if i-1>0 && griglia(i-1,j)==-1
                        griglia(i-1,j)=k;
                    end
                    % ne
                    if i-1>0 && j+1<m+1 && griglia(i-1,j+1)==-1
                        griglia(i-1,j+1)=k;
                    end
                    % e
                    if j+1<m+1 && griglia(i,j+1)==-1
                        griglia(i,j+1)=k;
                    end
                    % se
                    if i+1<n+1 && j+1<m+1 && griglia(i+1,j+1)==-1
                        griglia(i+1,j+1)=k;
                    end
                    % s
                    if i+1<n+1 && griglia(i+1,j)==-1
                        griglia(i+1,j)=k;
                    end
                    % sw
                    if i+1<n+1 && j-1>0 && griglia(i+1,j-1)==-1
                        griglia(i+1,j-1)=k;
                    end
                    % w
                    if j-1>0 && griglia(i,j-1)==-1
                        griglia(i,j-1)=k;
                    end
                end
            end
        end
    end

    % Rimozione degli ostacoli
    for o=5:4:length(ostacoli)
        for i=1:n
            for j=1:m
                if j>=ostacoli(o,1) && j<=ostacoli(o+1,1) && ...
                    i>=ostacoli(o+1,2) && i<=ostacoli(o+2,2)
                    griglia(j,i)=inf;
                end
            end
        end
    end

    % Rimozione dei lati della stanza
    griglia(1,:)=inf;griglia(:,1)=inf;
    griglia(n,:)=inf;griglia(:,m)=inf;

    % Plot
    [X,Y] = meshgrid(linspace(0,n,m),linspace(0,n,m));
    figure(); mesh(Y,X,griglia); axis equal;
    title('Approssimazione discreta del potenziale');
    xlabel('x');ylabel('y');zlabel('z');
    
    %% CALCOLO DEL PERCORSO
    path = start;
    PC = start;
    for i=1:n*m
        if PC==goal
            figure(2); hold on; title('Percorso con Discrete Potential Fields');
            x = path(:, 1); y = path(:, 2);
            plot(x, y, '-');
            return
        end
        
        [minValue, minX, minY] = findMinValueInSquare(griglia,PC(1),PC(2),3);
        
        if minValue<griglia(PC(1),PC(2))
            path=[path;minX minY];
        end
        
        if minValue==griglia(PC(1),PC(2))
            % per ogni cella ad un passo da minX,minY devo vedere quella di
            % valore minValue che non e' ancora stata visitata
            
            % Definizione degli offset per le otto celle ad un passo
            offsets = [-1, -1;-1, 0;-1, 1;0, -1;0, 1;1, -1;1, 0;1, 1];

            % Scorrimento delle otto celle ad un passo
            for j = 1:size(offsets, 1)
                newRow = PC(1) + offsets(j, 1);
                newCol = PC(2) + offsets(j, 2);
                
                cellValue = griglia(newRow, newCol);
                
                if cellValue==minValue
                    giaPassato=false;
                    for ics=1:size(path,1)
                        if path(ics,:)==[newRow newCol]
                            giaPassato=true;
                        end
                    end
                    if giaPassato==false
                        path=[path;newRow newCol];
                        break
                    end
                end
                if j == size(offsets, 1)
                     % Se si arriva a questo punto vuol dire che non è 
                     % stato trovato un punto in cui andare quindi bisogna
                     % fare backtracking
            
                     % vediamo lungo quale asse ci stiamo muovendo
                     if path(end,1)==path(end-1,1)
                         % ci stiamo muovendo lungo y
                         pos = size(path,1);
                         prec=[0 0];
                         while pos-1>0 && path(pos,1)==path(pos-1,1)
                             prec=path(pos,:);
                             path(pos,:)=[];
                             pos=pos-1;
                         end
                         if path(pos,2)>prec(2)
                             path=[path; path(pos,1) path(pos,2)+1];
                         else
                             path=[path; path(pos,1) path(pos,2)-1];
                         end
                     elseif path(end,2)==path(end-1,2)
                         % ci stiamo muovendo lungo x
                         pos = size(path,1);
                         prec=[0 0];
                         while pos-1>0 && path(pos,2)==path(pos-1,2)
                             prec=path(pos,:);
                             path(pos,:)=[];
                             pos=pos-1;
                         end
                         if path(pos,1)>prec(1)
                             path=[path; path(pos,1)+1 path(pos,2)];
                         else
                             path=[path; path(pos,1)-1 path(pos,2)];
                         end
                     end
                end
            end
        end
        
        PC=path(end,:);
        
    end
    
    %% Plot
    % Questo plot avviene solo se la funzione termina senza aver trovato un
    % percorso che va dal punto di start al punto di goal
    figure(2); hold on; title('Percorso con Discrete Potential Fields');
    x = path(:, 1); y = path(:, 2);
    plot(x, y, '-');
end

