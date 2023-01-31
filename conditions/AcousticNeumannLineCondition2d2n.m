classdef AcousticNeumannLineCondition2d2n < LineCondition2d2n
    %   AcousticNeumannLineCondition2d2n imposes a displacement boundary condition
    %   on an acoustic fluid
    
    %   (see: FRANCK 2008)
    
    methods
        function obj = AcousticNeumannLineCondition2d2n(id, nodeArray)
            requiredPropertyNames = ["FREQUENCY","NUMBER_GAUSS_POINT"];
            
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
        
        function NeumannVector = computeRHSContribution(obj)
            u(1) = obj.nodeArray(1).getPropertyValue('NORMAL_DISPLACEMENT');
            u(2) = obj.nodeArray(2).getPropertyValue('NORMAL_DISPLACEMENT');
            
            % no computation required, if normal displacement is zero
            if all(~abs(u))
                NeumannVector = u';
                return
            end
            
            % frequency must be initialized
            omega = obj.getPropertyValue('FREQUENCY');
            if isnan(omega)
                msg = ['AcousticNeumannLineCondition2d2n: No frequency defined', ...
                    ' at condition ', num2str(obj.getId())];
                e = MException('MATLAB:bm_mfem:conditionPropertyMissing',msg);
                throw(e);
            end
            
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            l = computeLength(obj.nodeArray(1).getCoords, ...
                    obj.nodeArray(2).getCoords);
            temp=sparse(2,1);
            
            for i=1:p
                eta=g(i);
                N = computeShapeFunction(obj, eta);
                temp=temp+(w(i)*transpose(N)*(l/2));
            end
            
            NeumannVector=omega^2*temp.*(transpose(u));
        end
        
        function ids = getDofList(obj)
            ids = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
        end
        
    end
end

