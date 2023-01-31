classdef PointLoadCondition2d1n < Condition
    %   POINTLOADCONDITION2D1N A point load condition for 2d
    %   The condition uses the nodal property POINT_LOAD
    
    methods
        function obj = PointLoadCondition2d1n(id, nodeArray)
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
            obj.dofNames = ["DISPLACEMENT_X", "DISPLACEMENT_Y"]; %do we need that?
        end
        
        function contrib = computeRHSContribution(obj)
            %computeRHSContribution Computes RHS contribution derived from
            %   the nodal property POINT_LOAD
            contrib = obj.nodeArray.getPropertyValue('POINT_LOAD');
            contrib = contrib(1:2);
            contrib = reshape(contrib,2,1);
        end
        
        function ids = getDofList(obj)
            ids = [obj.nodeArray.getDof('DISPLACEMENT_X'), ...
                obj.nodeArray.getDof('DISPLACEMENT_Y')];
        end
        
        function initialize(obj)

        end
        
    end
end

