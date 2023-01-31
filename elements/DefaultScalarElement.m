classdef DefaultScalarElement < Element
    %DEFAULTSCALARELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function initialize(obj); end
        function update(obj); end
        function barycenter(obj); end
        
        function dof = getDofList(obj)
            dof = Dof(Node(),[],[]);
            dof.setId(1);
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,length(obj.nodeArray.getDofArray()));
        end        
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            stiffnessMatrix = 0;
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            massMatrix = 0;
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            dampingMatrix = 0;
        end
    end
end

