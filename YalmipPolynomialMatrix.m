classdef YalmipPolynomialMatrix
    
    properties
        
        minDeg, 
        maxDeg, 
        varList, 
        coeffList, 
        monos, 
        polyRep,
        row,
        col, 
        
    end
    
    methods
        
        function obj = YalmipPolynomialMatrix(varListDef, maxDegreeDef, minDegreeDef, rowDef, colDef)
            
            obj.maxDeg = maxDegreeDef; 
            obj.minDeg = minDegreeDef;
            obj.varList = varListDef; 
            
            obj.polyRep = sdpvar(rowDef, colDef); 
            obj.row = rowDef; 
            obj.col = colDef; 
            obj.coeffList = [];
            obj.monos = []; 
            
            for i = 1:rowDef
                for j = 1:i
                    
                    [polyRepTemp, coeffListTemp, monosTemp] = polynomial(varListDef, maxDegreeDef, minDegreeDef); 
                    
                    obj.polyRep(i, j) = polyRepTemp; 
                    obj.polyRep(j, i) = polyRepTemp; 
                    obj.coeffList = [obj.coeffList; coeffListTemp]; 
                    obj.monos = [obj.monos; monosTemp]; 
                    
                end
            end
            
        end
        
        function matrixDot = lieDerivative(obj, f)
            
            matrixDot = sdpvar(obj.row, obj.col); 
            
            for i = 1:obj.row
                for j = 1:obj.col
                    matrixDot(i, j) = jacobian(obj.polyRep(i,j), obj.varList) * f; 
                end
            end
            
        end
        
    end
    
end