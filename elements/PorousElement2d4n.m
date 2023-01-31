classdef PorousElement2d4n < Element
    %   PorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        lengthX
        lengthY
    end
    
    methods
        % constructor
        function obj = PorousElement2d4n(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            obj@Element(super_args{:});
        end
        
        function c = barycenter(obj)
            %POROUSELEMENT3D4N.BARYCENTER returns the barycenter of the
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
            %POROUSELEMENT3D4N.CHECKCONVEXITY returns true, if the
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
        
        % member functions
        function [N_mat, N, Be, B, J] = computeShapeFunction(obj,xi,eta)
            % Shape Function and Derivatives
            N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];
            
            N_mat = sparse(2,8);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);
            
            % Coordinates of the nodes forming one element
            ele_coords = zeros(4,2);
            for i=1:4
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            % Jacobian
            J = N_Diff_Par * ele_coords;
            
            % Calculation of B-Matrix
            B=J\N_Diff_Par;
            Bx=B(1,1:4);
            By=B(2,1:4);
            
            Be=[Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
                0,By(1),0,By(2),0,By(3),0,By(4);
                By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
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
            x = [obj.nodeArray(1).getX + scaling * imag(obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_X', step)), ...
                obj.nodeArray(2).getX + scaling * imag(obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_X', step)), ...
                obj.nodeArray(3).getX + scaling * imag(obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_X', step)), ...
                obj.nodeArray(4).getX + scaling * imag(obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_X', step)),...
                obj.nodeArray(1).getX + scaling * imag(obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_X', step))];
            
            y = [obj.nodeArray(1).getY + scaling * imag(obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Y', step)), ...
                obj.nodeArray(2).getY + scaling * imag(obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_Y', step)), ...
                obj.nodeArray(3).getY + scaling * imag(obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_Y', step)), ...
                obj.nodeArray(4).getY + scaling * imag(obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_Y', step)),...
                obj.nodeArray(1).getY + scaling * imag(obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Y', step))];
            
            pl = line(x,y);
        end
        
        function update(obj)  
        end
        
        function [normalVector,dirX] = getNormalVector(obj,nodes)
            % getNormalVector returns normal vector on line defined by two
            % nodes 
            
            node1 = nodes(1).getCoords;
            node2 = nodes(2).getCoords;
            
            dirX = node2(1:2) - node1(1:2);
            dirX = dirX ./ norm(dirX);
            
            dirY = [-dirX(2) dirX(1)];
            
            tMat = [dirX; dirY];
            
            normalVector = tMat(2,:);
        end
        
    end
    
    % Defines element type for GIDOutput
    methods (Static)
        function o = getElementType()
            o = 'Quadrilateral';
        end
    end
    
end


