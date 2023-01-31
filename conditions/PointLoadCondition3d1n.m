classdef PointLoadCondition3d1n < Condition
    %   POINTLOADCONDITION3D1N A point load condition
    %   The condition uses the nodal property POINT_LOAD
    
    methods
        function obj = PointLoadCondition3d1n(id, nodeArray)
            %POINTLOADCONDITION3D1N Construct an instance of PointLoadCondition3d1n
            %   Parameters: id, nodeArray
            requiredPropertyNames = [];
            
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 1 && isa(nodeArray,'Node'))
                    error('problem with the nodes in condition %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            obj@Condition(super_args{:});
            obj.dofNames = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"]; %do we need that?
        end
        
        function contrib = computeRHSContribution(obj)
            %computeRHSContribution Computes RHS contribution derived from
            %   the nodal property POINT_LOAD
            contrib = obj.nodeArray.getPropertyValue('POINT_LOAD');
            contrib = reshape(contrib,3,1);
        end
        
        function ids = getDofList(obj)
            ids = [obj.nodeArray.getDof('DISPLACEMENT_X'), ...
                obj.nodeArray.getDof('DISPLACEMENT_Y'), ...
                obj.nodeArray.getDof('DISPLACEMENT_Z')];
        end
        
        function initialize(obj)

        end
        
    end
end

