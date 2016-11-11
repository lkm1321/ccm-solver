function matrixDot = getMatrixDot(matrix, f, vars)

    [row ,col] = size(matrix); 
    matrixDot = NaN(row, col); 
    
    for i = 1:row
        for j = 1:col
            matrixDot(i, j) = jacobian(matrix, vars) * f; 
        end
    end

end