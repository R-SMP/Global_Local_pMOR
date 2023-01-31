classdef FemModelPart < handle
    %FEMMODELPART Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        nodeArray
        elementArray
        conditionArray
        dofArray
        fixedDofs
        freeDofs
        prescribedDofs
        mProperties
        parentFemModel
    end
    
    methods
        
        function obj = FemModelPart(name, nodes, elements, conditions, properties, parentFemModel)
            obj.name = name;
            obj.nodeArray = nodes;
            obj.elementArray = elements;
            obj.conditionArray = conditions;
            if isa(properties,'PropertyContainer')
                obj.mProperties = properties;
            else
                obj.mProperties = PropertyContainer();
            end
            if nargin > 5
                obj.parentFemModel = parentFemModel;
            end
        end
        
        function initialize(obj)
            if ~isempty(obj.nodeArray)
                tmp = arrayfun(@(node) node.getDofArray, obj.nodeArray, 'UniformOutput', false);
                obj.dofArray = [tmp{:}];

                fixed = obj.dofArray.isFixed();
                obj.fixedDofs = obj.dofArray(fixed);
                obj.freeDofs = obj.dofArray(~fixed);
                
                prescribed = obj.dofArray.isPrescribed();
                obj.prescribedDofs = obj.dofArray(prescribed);
            end
            
            for e = obj.elementArray
                e.setProperties(obj.mProperties);
            end
            
            % (1a) assign properties to conditions
            for c = obj.conditionArray
                c.setProperties(obj.mProperties);
            end
        end
        
        function n = getName(obj)
            n = obj.name;
        end
        
        function n = getNodes(obj)
            n = obj.nodeArray;
        end
        
        function n = getNodeById(obj, id)
            n = Node.empty;
            if ~ isempty(obj.nodeArray)
                nodeIds = obj.nodeArray.getId();
                n = obj.nodeArray(id==nodeIds);
            end
            if isempty(n)
                msg = ['FemModelPart: Node with id ', num2str(id), ...
                    ' not found in modelpart ', obj.name];
                e = MException('MATLAB:bm_mfem:nodeNotFound',msg);
                throw(e);
            end
        end
        
        function e = getElements(obj)
            e = obj.elementArray;
        end
        
        function e = getElementById(obj, id)
            e = Element.empty;
            if ~ isempty(obj.elementArray)
                elementIds = obj.elementArray.getId();
                e = obj.elementArray(id==elementIds);
            end
            if isempty(e)
                msg = ['FemModelPart: Element with id ', num2str(id), ...
                    ' not found in modelpart ', obj.name];
                err = MException('MATLAB:bm_mfem:elementNotFound',msg);
                throw(err);
            end
        end
        
        function c = getConditions(obj)
            c = obj.conditionArray;
        end
        
        function c = getConditionById(obj, id)
            c = Condition.empty;
            if ~ isempty(obj.conditionArray)
                conditionIds = obj.conditionArray.getId();
                c = obj.conditionArray(id==conditionIds);
            end
            if isempty(c)
                msg = ['FemModelPart: Condition with id ', num2str(id), ...
                    ' not found in modelpart ', obj.name];
                err = MException('MATLAB:bm_mfem:conditionNotFound',msg);
                throw(err);
            end
        end
        
        function d = getDofArray(obj)
            d = obj.dofArray;
        end
        
        function [free, fix, prescribed] = getDofConstraints(obj)
            free = obj.freeDofs;
            fix = obj.fixedDofs;
            prescribed = obj.prescribedDofs;
        end
        
        function p = getParentFemModel(obj)
            p = obj.parentFemModel;
        end
        
        function mProperties = getProperties(femModel)
            mProperties = femModel.mProperties;
        end
        
        function setPropertyValue(obj, valueName, value)
            obj.mProperties.setValue(valueName, value);
        end
        
        function addNode(obj, node)
            if isa(node,'Node')
                obj.nodeArray = [obj.nodeArray node];
            else
                obj.nodeArray = [obj.nodeArray obj.parentFemModel.getNode(node)];
            end
        end
        
        function addElement(obj, element)
            if isa(element,'Element')
                obj.elementArray = [obj.elementArray element];
            else
                obj.elementArray = [obj.elementArray obj.parentFemModel.getElement(element)];
            end
        end
        
        function n = addNewNode(obj, id, x, y, z)
            if ~isempty(obj.parentFemModel)
                n = obj.parentFemModel.addNewNode(id, x, y, z);
                obj.nodeArray = [obj.nodeArray n];
            else
                msg = 'FemModelPart: No parent model defined';
                e = MException('MATLAB:bm_mfem:noParentModel',msg);
                throw(e);
            end
        end
        
        function e = addNewElement(obj, elementName, id, nodes, props)
            if ~isempty(obj.parentFemModel)
                if nargin == 4
                    e = obj.parentFemModel.addNewElement(elementName, id, nodes);
                elseif nargin == 5
                    e = obj.parentFemModel.addNewElement(elementName, id, nodes, props);
                else
                    msg = 'FemModelPart: Wrong number of input arguments';
                    err = MException('MATLAB:bm_mfem:invalidArguments',msg);
                    throw(err);
                end
                obj.elementArray = [obj.elementArray e];
            else
                msg = 'FemModelPart: No parent model defined';
                err = MException('MATLAB:bm_mfem:noParentModel',msg);
                throw(err);
            end
        end
        
        function e = addNewCondition(obj, conditionName, id, nodes, props)
            if ~isempty(obj.parentFemModel)
                if nargin == 4
                    e = obj.parentFemModel.addNewCondition(conditionName, id, nodes);
                elseif nargin == 5
                    e = obj.parentFemModel.addNewCondition(conditionName, id, nodes, props);
                else
                    msg = 'FemModelPart: Wrong number of input arguments';
                    err = MException('MATLAB:bm_mfem:invalidArguments',msg);
                    throw(err);
                end
                obj.conditionArray = [obj.conditionArray e];
            else
                msg = 'FemModelPart: No parent model defined';
                err = MException('MATLAB:bm_mfem:noParentModel',msg);
                throw(err);
            end
        end
    end
    
end

