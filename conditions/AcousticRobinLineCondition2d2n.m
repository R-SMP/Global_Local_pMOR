classdef AcousticRobinLineCondition2d2n < LineCondition2d2n
    %   AcousticRobinLineCondition2d2n imposes an impedance boundary condition
    %   on an acoustic fluid
    
    %   (see: FRANCK 2008)
    
    
    methods
        function obj = AcousticRobinLineCondition2d2n(id, nodeArray)
            requiredPropertyNames = ["NUMBER_GAUSS_POINT"];
            
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
            obj.dofNames = ('ACOUSTIC_PRESSURE');
        end
        
        function initialize(obj)
            % Default setting for NUMBER_GAUSS_POINT
            obj.props.setValue('NUMBER_GAUSS_POINT',2);
        end
        
        function admittanceMatrix = computeCContribution(obj)
            %computeCContribution Computes dampingMatrix contribution derived from
            %the Condition properties 'DENSITY_FLUID' and 'ADMITTANCE'
            a(1) = obj.nodeArray(1).getPropertyValue('ADMITTANCE');
            a(2) = obj.nodeArray(2).getPropertyValue('ADMITTANCE');
            A = (a(1)+a(2))/2;
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g] = returnGaussPoint(p);
            l = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(2).getCoords);
            admittanceMatrix = sparse(2,2);
            
            for i=1:p
                eta = g(i);
                N = computeShapeFunction(obj, eta);
                admittanceMatrix = admittanceMatrix+(w(i)*A*transpose(N)*N*(l/2));
            end
        end
        
        function ids = getDofList(obj)
            ids = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
        end
        
    end
end

