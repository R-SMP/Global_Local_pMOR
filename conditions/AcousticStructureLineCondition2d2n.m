classdef AcousticStructureLineCondition2d2n < LineCondition2d2n
    %   AcousticStructureLineCondition2d2n couples a fluid and a solid at a
    %   given coupling-boundary 
    
    %   DIRECTION == 1 if normal vector is pointing outward from acoustic
    %   domain (default)
    
    %   (see: RUMPLER 2012)
    
    methods
        function obj = AcousticStructureLineCondition2d2n(id, nodeArray)
            requiredPropertyNames = ["DIRECTION","NUMBER_GAUSS_POINT"];
            
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
                "DISPLACEMENT_Y"]);
        end
        
        function initialize(obj)
            % Default setting for DIRECTION of normal vector and NUMBER_GAUSS_POINT
            obj.props.setValue('DIRECTION',1);
            obj.props.setValue('NUMBER_GAUSS_POINT',2);
        end
        
        function couplingKMatrix = computeKContribution(obj)
            d = obj.getPropertyValue('DIRECTION');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g] = returnGaussPoint(p);
            tmat = getTransformationMatrix(obj);
            normalvector = tmat(2,:);
            temp = sparse(2,4);
            
            for i=1:p
                eta = g(i);
                [N, N_mat, J] = computeShapeFunction(obj, eta);
                temp = temp+(w(i)*transpose(N)*normalvector*N_mat*J);
            end
            
            couplingKMatrix = sparse(6,6);
            couplingKMatrix(1:4,5:6) = d * (-1) * temp.';
        end
        
        function couplingMMatrix = computeMContribution(obj)
            d = obj.getPropertyValue('DIRECTION');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g] = returnGaussPoint(p);
            tmat = getTransformationMatrix(obj);
            normalvector = tmat(2,:);
            temp = sparse(2,4);
            
            for i=1:p
                eta = g(i);
                [N, N_mat, J] = computeShapeFunction(obj, eta);
                temp = temp+(w(i)*transpose(N)*normalvector*N_mat*J);
            end
            
            couplingMMatrix = sparse(6,6);
            couplingMMatrix(5:6,1:4) = d * temp;
        end
        
        function ids = getDofList(obj)
            ids([1 3]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            ids([2 4]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            ids([5 6]) = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
        end
        
    end
end

