classdef LineCondition2d2n < Condition
    %   LineCondition2d2n is an abstract class to model a condition that
    %   applies on a linear boundary or coupling interface
    
     properties (Access = protected)
        length
     end
    
    methods
        function obj = LineCondition2d2n(id, nodeArray, requiredPropertyNames)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            obj@Condition(super_args{:});
        end
        
        function [N, N_mat, J] = computeShapeFunction(obj, eta)
            % Computes ShapeFunction at input Gaußpoint
            N = [(1-eta)/2    (1 + eta)/2];
            N_mat = [ N(1), 0, N(2), 0; 0, N(1), 0, N(2)];
            
            [c1, c2] = getTransformedCoords(obj);
            
            eleCoords = [c1(1); c2(1)];
            J = [-1/2 1/2] * eleCoords; 
        end
        
        function len = getLength(obj)
            % Computes length of an Element - for choosen shape functions the JacobiDeterminant equals half the element length 
            if obj.length == 0
                obj.length = computeLength(obj.nodeArray(1).getCoords, ...
                    obj.nodeArray(2).getCoords);
            end
            len = obj.length;
        end
        
        function tMat = getTransformationMatrix(obj)
            %Computes the TransferMatrix that transforms the global
            %coordinates of the boundary into local coordinates
            node1 = obj.nodeArray(1).getCoords;
            node2 = obj.nodeArray(2).getCoords;
            
            dirX = node2(1:2) - node1(1:2);
            dirX = dirX ./ norm(dirX);

            dirY = [-dirX(2) dirX(1)];

            tMat = [dirX; dirY];
        end
        
        function [c1, c2] = getTransformedCoords(obj)
            %Computes the local Coordinates of the boundary
            tMat = getTransformationMatrix(obj);
            node1 = obj.nodeArray(1).getCoords;
            node2 = obj.nodeArray(2).getCoords;
            
            c1 = tMat * (node1(1:2))';
            c2 = tMat * (node2(1:2))';
        end
        
    end
end

