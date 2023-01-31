classdef AcousticElement2d2n < Element
    %   AcousticElement2d2n Summary of this class goes here
    %   Detailed explanation goes here
   
    %   Can be used if acoustic pressure excitation should be imposed on 
    %   element that does not contain acoustic pressure as DOF
    properties (Access = protected)
        lengthX
        lengthY
        lenghtZ
    end
    
    methods
        % constructor
        function obj = AcousticElement2d2n(id, nodeArray)
            requiredProperties = cellstr(["YOUNGS_MODULUS","DENSITY","NUMBER_GAUSS_POINT"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredProperties};
            end
            
            % call the super class constructor
            obj@Element(super_args{:});
            obj.dofNames = cellstr('ACOUSTIC_PRESSURE');
        end
        
        
        %Initialization
        function initialize(obj)
            % Default setting for parameters of fluid (AIR) and
            % NUMBER_GAUSS_POINT
            obj.eProperties.setValue('DENSITY',1.21);
            obj.eProperties.setValue('YOUNGS_MODULUS',142771);
            obj.eProperties.setValue('NUMBER_GAUSS_POINT',2);
        end
        
        % member functions
        function [N, B, J] = computeShapeFunction(obj,eta)
            % Shape Function and Derivatives
            N = [(1-eta)/2    (1+eta)/2];
            N_Diff_Par = [-1/2    1/2];

            [c1, c2] = getTransformedCoords(obj);
            
            eleCoords = [c1(1); c2(1)];
            J = [-1/2 1/2] * eleCoords; 
            
            % Calculation of B-Matrix
            B=J\N_Diff_Par;
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
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            rho = obj.tryGetPropertyValue('DENSITY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            stiffnessMatrix=sparse(2,2);
            
            
            for i=1:p
                eta=g(i);
                [~, B, J] = computeShapeFunction(obj, eta);
                stiffnessMatrix=stiffnessMatrix+(w(i)*(1/rho)*det(J)*transpose(B)*B);
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            E = obj.tryGetPropertyValue('YOUNGS_MODULUS');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            massMatrix=sparse(2,2);
            
            
            for i=1:p
                eta=g(i);
                [N, ~, J] = computeShapeFunction(obj, eta);
                massMatrix=massMatrix+(w(i)*(1/E)*transpose(N)*N*det(J));
            end
            
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            dampingMatrix=zeros(2,2);
        end

        function dofs = getDofList(element)
            dofs([1 2]) = element.nodeArray.getDof('ACOUSTIC_PRESSURE');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,2);
            vals([1 2]) = element.nodeArray.getDofValue('ACOUSTIC_PRESSURE',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,2);
            [~, vals([1 2]), ~] = element.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,2);
            [~, ~, vals([1 2])] = element.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function update(obj)
        end
        
        function barycenter(obj)
        end
        
    end
    
    methods (Static)
        function o = getElementType()
            o = 'Linear';
        end
    end
    
end


