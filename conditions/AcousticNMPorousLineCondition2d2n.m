classdef AcousticNMPorousLineCondition2d2n < LineCondition2d2n
    %   AcousticNMPorousLineCondition2d2n couples an acoustic fluid and a porous medium (new mixed formulation) at a
    %   given coupling-boundary
    
    %   DIRECTION == 1 if normal vector is pointing outward from porous
    %   domain (default)
    
    %   (see: ATALLA 2001)
    
    methods
        function obj = AcousticNMPorousLineCondition2d2n(id, nodeArray)
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
            couplingKMatrix(3:6,1:2) = (d) * (temp.');
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
            couplingMMatrix(1:2,3:6) = -d * temp;
        end
        
        function ids = getDofList(obj)
            ids([1 2]) = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
            ids([3 5]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            ids([4 6]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
    end
end

