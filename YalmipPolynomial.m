classdef YalmipPolynomial
    
    properties
        
        minDeg, 
        maxDeg, 
        varList, 
        coeffList, 
        monos, 
        polyRep,
        
    end
    
    methods
        function obj = YalmipPolynomial(varListDef, maxDegreeDef, minDegreeDef)
            
            obj.maxDeg = maxDegreeDef; 
            obj.minDeg = minDegreeDef;
            obj.varList = varListDef; 
            [obj.polyRep, obj.coeffList, obj.monos] = polynomial(varListDef, maxDegreeDef, minDegreeDef); 
            
        end
    end
    
end