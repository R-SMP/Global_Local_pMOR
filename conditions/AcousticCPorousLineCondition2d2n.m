classdef AcousticCPorousLineCondition2d2n < LineCondition2d2n
    %   AcousticCPorousLineCondition2d2n couples a fluid and a porous medium (classical formulation) at a
    %   given coupling-boundary 
    
    %   DIRECTION == 1 if normal vector is pointing outward from acoustic
    %   domain (default)
    
    %   (see: RUMPLER 2012)
    
    methods
        function obj = AcousticCPorousLineCondition2d2n(id, nodeArray)
            requiredPropertyNames = ["POROSITY","DIRECTION","NUMBER_GAUSS_POINT"];
            
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    error('problem with the nodes in condition %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            obj@LineCondition2d2n(super_args{:});
            obj.dofNames = cellstr(["ACOUSTIC_PRESSURE", "DISPLACEMENT_X", ...
                "DISPLACEMENT_Y", "DISPLACEMENT_FLUID_X", "DISPLACEMENT_FLUID_Y"]);
        end
        
        function initialize(obj)
            % Default setting for DIRECTION of normal vector and NUMBER_GAUSS_POINT
            obj.props.setValue('DIRECTION',1);
            obj.props.setValue('NUMBER_GAUSS_POINT',2);
        end
        
        function CouplingKMatrix = computeKContribution(obj)
            d = obj.getPropertyValue('DIRECTION');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            POROSITY = obj.getPropertyValue('POROSITY');
            [w,g]=returnGaussPoint(p);
            tmat=getTransformationMatrix(obj);
            normalvector = tmat(2,:);
            temp=sparse(2,4);
            
            for i=1:p
                eta=g(i);
                [N, N_mat, J] = computeShapeFunction(obj, eta);
                temp=temp+(w(i)*transpose(N)*normalvector*N_mat*J);
            end
            
            CouplingKMatrix = sparse(10,10);
            CouplingKMatrix(3:6,1:2) = d * (-1) * (1 - POROSITY) * temp.';
            CouplingKMatrix(7:10,1:2) = d * (-1) * POROSITY * temp.';
        end
        
        function CouplingMMatrix = computeMContribution(obj)
            d = obj.getPropertyValue('DIRECTION');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            POROSITY = obj.getPropertyValue('POROSITY');
            tmat=getTransformationMatrix(obj);
            normalvector = tmat(2,:);
            temp=sparse(2,4);
            
            for i=1:p
                eta=g(i);
                [N, N_mat, J] = computeShapeFunction(obj, eta);
                temp=temp+(w(i)*transpose(N)*normalvector*N_mat*J);
            end
            
            CouplingMMatrix = sparse(10,10);
            CouplingMMatrix(1:2,3:6) = d * (1 - POROSITY) * temp;
            CouplingMMatrix(1:2,7:10) = d * POROSITY * temp;
        end
        
        function ids = getDofList(obj)
            ids([1 2]) = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
            ids([3 5]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            ids([4 6]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            ids([7 9]) = obj.nodeArray.getDof('DISPLACEMENT_FLUID_X');
            ids([8 10]) = obj.nodeArray.getDof('DISPLACEMENT_FLUID_Y');
        end
        
    end
end

