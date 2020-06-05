classdef PolyFit

    properties
        crdx
        crdy
        data
        consCoef % all coefficients that are fully determioned via the user-input constraints.
        freeCoef
        polyDegree
    end

    methods

        function self = PolyFit(polyDegree, crdx, crdy, data)

            crdxLen = length(crdx);
            crdyLen = length(crdy);
            if crdxLen ~= crdyLen
                error( "The length of crdx and crdy must be equal. These two represent the boundaries of the polynomial fit, " ...
                     + "or more accurately, the points through which the polynomial fit must pass." ...
                     );
            end

            if polyDegree < crdxLen
                error( "degree cannot be less than the length of crdx or crdy. This would result in an over-deterministic system of equations." );
            end

        end
        
        function polyFunc = getPolyFunc(self,freeCoef)

            consCoef = zeros(crdxLen,1);
            for i = 1:crdxLen
                
                consCoef(i) = crdy(i) - sum( freeCoef(i: .* crdx.^(i) );
                \Sum_{i=1}^{i=N} a_i * x0^i;
                a0 = y0 - \Sum_{i=1}^{i=N} a_i * x0^i
                
            end
            
        end
    
    end
    
end