function [yalmipMatrix, matrixDot] = getPolynomialMatrix(variables, maxDeg, minDeg, row, col)

    yalmipMatrix = NaN(row, col); 
    
    for i = 1:row
        for j = 1:col
                yalmipMatrix(i, j) = polynomial(variables, maxDeg, minDeg); 
        end
    end

end 