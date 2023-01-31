classdef Condition < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %CONDITION Base class for conditions
    
    properties (Access = private)
        id
    end
    
    properties (Access = protected)
        nodeArray
        dofNames %do we need that?
        props
        requiredPropertyNames
    end
    
    methods
        function obj = Condition(id, nodeArray, requiredPropertyNames)
            %CONDITION Construct an instance of Condition
            %   Parameters: id, nodeArray, requiredPropertyNames
            
            obj.props = PropertyContainer();
            
            if nargin == 0
                obj.id = -1;
            elseif nargin == 3
                obj.id = id;
                obj.nodeArray = nodeArray;
                
                for ii = 1:length(requiredPropertyNames)
                    obj.props.addValue(requiredPropertyNames{ii},NaN);
                end
                
                obj.initialize();
                
            else
                error('wrong number of input arguments')
            end
            
        end
        
    end
    
    methods (Sealed)
        
        function id = getId(obj)
            id = zeros;
            for ii = 1:length(obj)
                id(ii) = obj(ii).id;
            end
        end
        
        function n = getNodes(obj)
            n = obj.nodeArray;
        end
        
        function prop = getProperties(obj)
            prop = obj.props;
        end
        
        function value = getPropertyValue(obj, valueName)
           value = obj.props.getValue(valueName); 
        end
        
        function setProperties(obj, newProps)
            for ii = 1:length(obj)
                names = newProps.getValueNames();
                for jj = 1:length(names)
                    obj(ii).props.setValue(names{jj}, newProps.getValue(names{jj}));
                end
            end
        end
        
        function setPropertyValue(obj, valueName, value)
            for ii = 1:length(obj)
                obj(ii).props.setValue(valueName, value);
            end
        end
    end
    
    methods (Abstract)
        getDofList(obj)
        initialize(obj)
    end
    
    methods
        
        function c = computeRHSContribution(obj)
            n = length(obj.getDofList());
            c = sparse(n,1);
        end
        
        function c = computeKContribution(obj)
            n = length(obj.getDofList());
            c = sparse(n,n);
        end
        
        function c = computeCContribution(obj)
            n = length(obj.getDofList());
            c = sparse(n,n);
        end
        
        function c = computeMContribution(obj)
            n = length(obj.getDofList());
            c = sparse(n,n);
        end
        
    end
    
    methods (Static, Sealed, Access = protected)
        function obj = getDefaultScalarElement
            obj = DefaultScalarCondition;
        end
    end
end

