function [minValue, minX, minY] = findMinValueInSquare(matrix, x, y, s)
%FINDMINVALUEINSQUARE
%   Restituisce le coordinate del valore minimo tra quelli presenti nel 
%   quadrato che circonda il punto (x, y) della matrice matrix
%   matrix = matrice nxm
%   x, y = coordinate del punto
%   s = dimensione del quadrato circoscritto

    [rows, cols] = size(matrix);
    
    % Calcola gli indici di inizio e fine per il quadrato circoscritto
    startRow = max(1, x - floor(s/2));
    endRow = min(rows, x + ceil(s/2) - 1);
    startCol = max(1, y - floor(s/2));
    endCol = min(cols, y + ceil(s/2) - 1);
    
    % Inizializza il valore minimo con il primo elemento nel quadrato
    minValue = matrix(startRow, startCol);
    minX = startRow;
    minY = startCol;
    
    % Cerca il valore minimo all'interno del quadrato
    for i = startRow:endRow
        for j = startCol:endCol
            if matrix(i, j) < minValue
                minValue = matrix(i, j);
                minX = i;
                minY = j;
            end
        end
    end
    
end

