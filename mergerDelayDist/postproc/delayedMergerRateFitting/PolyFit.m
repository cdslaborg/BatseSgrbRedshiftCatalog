classdef PolyFit

    properties
        %data
        dofmax
        dof
        crdx
        crdy
        crdxLen
        crdyLen
        consCoef % all coefficients that are fully determioned via the user-input constraints.
        %freeCoef
        polyDegree
    end

    methods

        function self = PolyFit(polyDegree, crdx, crdy)

            %self.data = data;
            self.crdx = crdx;
            self.crdy = crdy;
            self.crdxLen = length(crdx); % number of constraints
            self.crdyLen = length(crdy); % number of constraints
            self.polyDegree = polyDegree;
            self.dofmax = polyDegree + 1;
            self.dof = self.dofmax - self.crdxLen;

            if self.crdxLen ~= self.crdyLen
                error( "The length of crdx and crdy must be equal. These two represent the boundaries of the polynomial fit, " ...
                     + "or more accurately, the points through which the polynomial fit must pass." ...
                     );
            end

            if self.polyDegree < self.crdxLen
                error( "degree cannot be less than the length of crdx or crdy. This would result in an over-deterministic system of equations." );
            end

        end
        
        function polyFunc = get(self,freeCoef,x)
            % the length of freeCoef is ( 1 : self.polyDegree - self.crdxLen )

            if length(freeCoef) ~= self.dof
                error("length(freeCoef) ~= self.dof");
            end
            
            % get the fixed coefficients resulting from the constraints

            self.consCoef = zeros(1,self.crdxLen);
            for k = 1:self.crdxLen

                if k==1
                    sumConsCoefTerm = 0;
                else
                    sumConsCoefTerm = sum( self.consCoef(1:k-1) .* self.crdx(k).^(0:k-2) );
                end
                
disp("k="+string(k));
disp( "freeCoef(1:self.dof)=" + string( freeCoef(1:self.dof) ) );
disp( "self.crdx(k).^(k+1:self.polyDegree)=" + string( self.crdx(k).^(k+1:self.polyDegree) ) );
                sumFreeCoefTerm = sum( freeCoef(self.dof-k+1:self.dof) .* self.crdx(k).^(k+1:self.polyDegree) );
sumFreeCoefTerm

                self.consCoef(k) = ( 1 / self.crdx(k) )^(k-1) ...
                                 * ( self.crdy(k) - sumFreeCoefTerm - sumConsCoefTerm );

            end

            % compute the polynomial

            vector1 = self.consCoef .* x.^(0:self.crdxLen-1);
            vector2 = freeCoef .* x.^(self.crdxLen:self.polyDegree-1);
            polyFunc = sum( vector1(:) ) + sum( vector2(:) );

        end
    
    end
    
end