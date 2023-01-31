classdef SeptumMPorousLineCondition2d2n < LineCondition2d2n
    %   SeptumMPorousLineCondition2d2n couples an imprevious septum and a porous medium (mixed formulation) at a
    %   given coupling-boundary;
    %   External forces are not taken into account in this coupling strategy
    %   but can be apllied as PointLoadCondition;
    %   SURFACE_DENSITY equals DENSITY*THICKNESS
    
    %   DIRECTION == 1 if normal vector is pointing outward from porous
    %   domain (default)
    
    %   (see: DEBERGUE 1999)
    
    methods
        function obj = SeptumMPorousLineCondition2d2n(id, nodeArray)
            requiredPropertyNames = ["SURFACE_DENSITY","DIRECTION","NUMBER_GAUSS_POINT"];
            
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
        
        function CouplingKMatrix = computeKContribution(obj)
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
            
            CouplingKMatrix = sparse(6,6);
            CouplingKMatrix(3:6,1:2) = d * (-temp.');
        end
        
        function CouplingMMatrix = computeMContribution(obj)
            d = obj.getPropertyValue('DIRECTION');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            tmat = getTransformationMatrix(obj);
            normalvector = tmat(2,:);
            m = obj.getPropertyValue('SURFACE_DENSITY');
            temp = sparse(2,4);
            temp1=sparse(4,4);
            
            for i=1:p
                eta=g(i);
                [N, N_mat, J] = computeShapeFunction(obj, eta);
                temp = temp+(w(i)*transpose(N)*normalvector*N_mat*J);
                temp1=temp1+(w(i)*m*transpose(N_mat)*N_mat*J);
            end
            
            CouplingMMatrix = sparse(6,6);
            CouplingMMatrix(1:2,3:6) = d * temp;
            CouplingMMatrix(3:6,3:6) = temp1;
        end
        
        function ids = getDofList(obj)
            ids([1 2]) = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
            ids([3 5]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            ids([4 6]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
    end
end

