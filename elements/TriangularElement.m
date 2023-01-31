classdef TriangularElement < Element
    %TRIANGULARELEMENT Base class for triangular elements
    %   Detailed explanation goes here
    
    properties (Access = protected)
    end
    
    methods
        % Constructor
        function obj = TriangularElement(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            
            obj@Element(super_args{:});
        end

        function c = barycenter(obj)
            c(1) = mean(obj.nodeArray.getX());
            c(2) = mean(obj.nodeArray.getY());
            c(3) = mean(obj.nodeArray.getZ());
        end
        
        function pl = draw(obj)   
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ... 
                 obj.nodeArray(3).getX, obj.nodeArray(1).getX];
             
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ... 
                 obj.nodeArray(3).getY, obj.nodeArray(1).getY];
             
            z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ... 
                 obj.nodeArray(3).getZ, obj.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end

        function update(obj)
        end
        
        function [sepTrig, sepCarth] = getStressEvaluationPoints(obj)
            nGaussPoints = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            sepTrig = zeros(nGaussPoints,3);
            for xi = 1 : nGaussPoints
                [~,g] = returnGaussPointTrig(nGaussPoints, xi);
                sepTrig(xi,:) = g;
            end
            
            if nGaussPoints == 1
                sepCarth = [1/3 1/3];
            elseif nGaussPoints == 3
                sepCarth = [1/6 1/6; 2/3 1/6; 1/6 2/3];
            else
                 msg = [class(obj), ': Carthesian coordinates for this ' ...
                     ' number of Gauss points not available.'];
                e = MException('MATLAB:bm_mfem:missingProperty',msg);
                throw(e);
            end
            
        end
        
        function computeElementStress(obj, step)
            %COMPUTEELEMENTSTRESS Compute the elemental stress tensor
            %   The stresses are evaluated at the inner Gauss points
            %   Only for two-dimensional triangles!
            if nargin==1; step = 'end'; end
            EModul = obj.getPropertyValue('YOUNGS_MODULUS');
            prxy = obj.getPropertyValue('POISSON_RATIO');
            sep = obj.getStressEvaluationPoints();
            s = zeros(3,3,size(sep,1));
            
            % Moment-Curvature Equations
            D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
            % Material Matrix D
            D = D * (EModul / (1-prxy^2));
            
            for ii = 1:size(sep,1)
                [~, ~, B, ~] = obj.computeShapeFunction(sep(ii,:));
                disp = obj.getValuesVector(step);
                strain = B * disp';
                stress = D * strain;
                s(1,1,ii) = stress(1);
                s(2,2,ii) = stress(2);
                s(1,2,ii) = stress(3);
                s(2,1,ii) = stress(3);
            end
            
            obj.setPropertyValue('STRESS_TENSOR', s);            
        end
        
    end
    
    methods (Static)
        function o = getElementType()
            o = 'Triangle';
        end
    end
    
end

