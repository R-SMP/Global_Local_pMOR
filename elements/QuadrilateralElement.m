classdef (Abstract) QuadrilateralElement < Element
    %QUADRILATERALELEMENT Base class for quadrilateral elements
    %   Detailed explanation goes here
    
    properties (Access = protected)
        lengthX
        lengthY
    end
    
    methods
        % Constructor
        function quadrilateralElement = QuadrilateralElement(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            
            quadrilateralElement@Element(super_args{:});
        end
        
        function c = barycenter(obj)
            %QUADRILATERALELEMENT.BARYCENTER returns the barycenter of the
            %element in the x-y plane
            diag1X = [obj.nodeArray(1).getX() obj.nodeArray(3).getX()];
            diag1Y = [obj.nodeArray(1).getY() obj.nodeArray(3).getY()];
            diag2X = [obj.nodeArray(2).getX() obj.nodeArray(4).getX()];
            diag2Y = [obj.nodeArray(2).getY() obj.nodeArray(4).getY()];
            
            try
                [c(1),c(2)] = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            catch e
                if (strcmp(e.identifier,'MATLAB:subsassigndimmismatch'))
                    msg = ['Barycenter of element ', ...
                        num2str(obj.getId()), ...
                        ' could not be computed. ', ...
                        'Check the element for convexity.'];
                    causeException = MException('MATLAB:bm_mfem:barycenter',msg);
                    e = addCause(e,causeException);
                elseif (strcmp(e.identifier,'MATLAB:UndefinedFunction'))
                    msg = ['Function ''polyxpoly'' not found. ', ...
                        'Please install the Mapping Toolbox.'];
                    causeException = MException('MATLAB:bm_mfem:toolboxnotfound',msg);
                    e = addCause(e,causeException);
                end
                rethrow(e);
            end
        end
        
        function c = checkConvexity(obj)
            %QUADRILATERALELEMENT.CHECKCONVEXITY returns true, if the
            %element is convex
            diag1X = [obj.nodeArray(1).getX() obj.nodeArray(3).getX()];
            diag1Y = [obj.nodeArray(1).getY() obj.nodeArray(3).getY()];
            diag2X = [obj.nodeArray(2).getX() obj.nodeArray(4).getX()];
            diag2Y = [obj.nodeArray(2).getY() obj.nodeArray(4).getY()];
            
            c = 0;
            try
                res = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            catch e
                if (strcmp(e.identifier,'MATLAB:subsassigndimmismatch'))
                    c = 0;
                elseif (strcmp(e.identifier,'MATLAB:UndefinedFunction'))
                    msg = ['Function ''polyxpoly'' not found. ', ...
                        'Please install the Mapping Toolbox.'];
                    causeException = MException('MATLAB:bm_mfem:toolboxnotfound',msg);
                    e = addCause(e,causeException);
                    rethrow(e);
                end
                
            end
            
            if ~isempty(res); c = 1; end
        end
        
        function pl = draw(obj)
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(4).getX,...
                obj.nodeArray(1).getX];
            
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(4).getY, ...
                obj.nodeArray(1).getY];
            
            if(all(obj.getNodes().getDimension == 3))
                z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                    obj.nodeArray(3).getZ, obj.nodeArray(4).getZ, ...
                    obj.nodeArray(1).getZ];
                
                pl = line(x,y,z);
            else
                pl = line(x,y);
            end
            
        end
        
        function pl = drawDeformed(obj, step, scaling)
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step),...
                obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
            
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            
            try
                z = [obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step),...
                    obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step)];
                
                pl = line(x,y,z);
            catch e
                if strcmp(e.identifier,'MATLAB:Containers:Map:NoKey')
                    pl = line(x,y);
                else
                    throw(e);
                end
            end
            
        end
        
        function update(obj)
        end
        
        function sep = getStressEvaluationPoints(obj, location)
            if nargin==1; location = 'onGaussPoints'; end
            if ~strcmp(location, 'onGaussPoints') && ~strcmp(location, 'onPoints')
                error('no')
            end
            
            if strcmp(location, 'onGaussPoints')
                nGaussPoints = obj.getPropertyValue('NUMBER_GAUSS_POINT');
                [~,g] = returnGaussPoint(nGaussPoints);
                sep = zeros(nGaussPoints^2,2);
                counter = 1;
                for xi = 1 : nGaussPoints
                    for eta = 1 : nGaussPoints
                        sep(counter,:) = [g(xi) g(eta)];
                        counter = counter + 1;
                    end
                end
            elseif strcmp(location, 'onPoints')
                sep = obj.stressPoints();
            end
            
        end
        
        function computeElementStress(obj, step)
            %COMPUTEELEMENTSTRESS Compute the elemental stress tensor
            %   The stresses are evaluated at the inner Gauss points
            %   Only for two-dimensional quads!
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
                [~, ~, B, ~] = obj.computeShapeFunction(sep(ii,1),sep(ii,2));
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
            o = 'Quadrilateral';
        end
    end
    
    methods (Access = protected)
        function [T, mid] = getTransformationMatrix(obj)
            %GETTRANSFORMATIONMATRIX returns the transformation matrix T of
            %the element with respect to its center point MID
            c1 = obj.nodeArray(1).getCoords();
            c2 = obj.nodeArray(2).getCoords();
            c3 = obj.nodeArray(3).getCoords();
            c4 = obj.nodeArray(4).getCoords();
            mid = 0.25 .* (c1 + c2 + c3 + c4);
            
            d13 = c3 - c1;
            d24 = c4 - c2;
            
            e3 = cross(d13,d24);
            e3 = e3/norm(e3);
            
            e1 = c2 - c1;
            e1 = e1 - dot(e1,e3)*e3;
            e1 = e1/norm(e1);
            
            e2 = cross(e3,e1);
            e2 = e2/norm(e2);
            
            T = sparse([e1;e2;e3]);
        end
    end
    
end
