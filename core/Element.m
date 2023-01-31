classdef (Abstract) Element < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %ELEMENT The element class
    %   Abstract base class for all element implementations
    %
    %   OBJ = Element( ID, NODEARRAY, REQUIREDPROPERTYNAMES)
    %       ID: unique id of the element.
    %       NODEARRAY: all nodes of the elements. Have to be class NODE
    %       REQUIREDPROPERTYNAMES: list of properties, which are required
    %           for a static analysis. These properties must have a value 
    %           which is not zero or NaN.
    
    properties (Access = private)
        id = -1
    end
    properties (Access = protected)
        nodeArray
        dofNames
        eProperties
        requiredPropertyNames
    end
    
    methods
        % constructor
        function obj = Element(id, nodeArray, requiredPropertyNames)
            
            if nargin == 0
                obj.id = -1;
            elseif nargin == 3
                obj.id = id;
                obj.eProperties = PropertyContainer();
                obj.nodeArray = nodeArray;
                obj.requiredPropertyNames = requiredPropertyNames;
                
                for ii = 1:length(requiredPropertyNames)
                    [~, default] = checkPropertyName(requiredPropertyNames{ii});
                    obj.eProperties.addValue(requiredPropertyNames{ii}, default);
                end
                
                obj.initialize();
            end
        end
        
    end
    
    methods (Abstract)
        initialize(element)
        update(element)     % update properties after e.g. nodes changed
        barycenter(element)
        computeLocalStiffnessMatrix(element)
%         computeLocalForceVector(element)
        getDofList(element)
        getValuesVector(element, step)
    end
    
    methods (Sealed)
        % getter functions
        function id = getId(element)
            id = zeros;
            for ii = 1:length(element)
                id(ii) = element(ii).id;
            end
        end
        
        function prop = getProperties(element)
            prop = element.eProperties;
        end
        
        function value = getPropertyValue(element, valueName)
           value = element.eProperties.getValue(valueName); 
        end
        
        function nodes = getNodes(elements)
            nodes = Node.empty;
            for ii = 1:length(elements)
                nodes(ii,:) = elements(ii).nodeArray;
            end
        end
        
%         function dofs = getDofs(element)
%            dofs = Dof.empty;
%            nodes = element.nodeArray;
%            for ii = 1:length(nodes)
%               dofs = [dofs nodes(ii).getDofArray];
%            end
%         end
        
        function setProperties(obj, props)
            for ii = 1:length(obj)
                names = props.getValueNames();
                for jj = 1:length(names)
                    if obj(ii).eProperties.hasValue(names{jj})
                        obj(ii).eProperties.setValue(names{jj}, props.getValue(names{jj}));
                    else
                        obj(ii).eProperties.addValue(names{jj}, props.getValue(names{jj}));
                    end
                end
            end
        end
        
        function setPropertyValue(elements, valueName, value)
            for ii = 1:length(elements)
                elements(ii).eProperties.setValue(valueName, value);
            end
        end
        
        function addProperty(elements, valueName, value)
            for ii = 1:length(elements)
                elements(ii).eProperties.addValue(valueName, value);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function cp = copyElement(obj)
           cp = copyElement@matlab.mixin.Copyable(obj);
           obj.id = obj.id + 100;
        end
        
%         function addDofs(element, dofNames)
%             for itNode = 1:length(element.nodeArray)
%                 nodalDofs(1, length(dofNames)) = Dof;
%                 for itDof = 1:length(dofNames)
%                     newDof = Dof(element.nodeArray(itNode),0.0,dofNames(itDof));
%                     nodalDofs(itDof) = newDof;
%                 end
%                 element.nodeArray(itNode).setDofArray(nodalDofs);
%             end
%         end
        
        function addDofsToSingleNode(element, node)
            nodalDofs(1, length(element.dofNames)) = Dof;
            for itDof = 1:length(element.dofNames)
                newDof = Dof(node,0.0,element.dofNames(itDof));
                nodalDofs(itDof) = newDof;
            end
            node.setDofArray(nodalDofs);
        end
        
    end
    
    methods
        
        function check(element)
        %CHECK checks, if all dofs and properties required by the element
        %   are available. If a property is missing, it will be initialized
        %   with 0.
            
            %check the dofs
            for iNode = 1:length(element.nodeArray)
                cNode = element.nodeArray(iNode);
                availableDofNames = arrayfun(@(dof) dof.getValueType, cNode.getDofArray, 'UniformOutput',false);
                diff = setxor(element.dofNames', availableDofNames);
                if ~ isempty(diff)
                    missingDofs = setdiff(element.dofNames', availableDofNames);
                    
                    if ~isempty(missingDofs)
                        msg = ['Element: The following dofs are missing', ...
                            ' at node ', num2str(cNode.getId), ': ', ...
                            char(strjoin(missingDofs,', '))];
                        e = MException('MATLAB:bm_mfem:elementalDofMissing',msg);
                        throw(e);
                    end
                    
                end
            end
            
            for ii = 1:length(element.requiredPropertyNames)
                prop = element.getPropertyValue(element.requiredPropertyNames{ii});
                if any(isnan(prop))
                    msg = [class(element), ': The required property \"', ...
                        element.requiredPropertyNames{ii}, '\" in element ', ...
                        num2str(element.getId), ' is not initialized. '];
                    e = MException('MATLAB:bm_mfem:requiredPropertyMissing',msg);
                    throw(e);
                end
            end

        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            %COMPUTELOCALDAMPINGMATRIX returns the damping matrix for
            %   proportional (Rayleigh) damping
            props = obj.getProperties;
            
            if props.hasValue('RAYLEIGH_ALPHA') && props.hasValue('RAYLEIGH_BETA')
                alpha = props.getValue('RAYLEIGH_ALPHA');
                beta = props.getValue('RAYLEIGH_BETA');
                dampingMatrix = alpha * obj.computeLocalMassMatrix + ...
                    beta * obj.computeLocalStiffnessMatrix;
            else
                ndofs = length(obj.getDofList());
                dampingMatrix = sparse(ndofs, ndofs);
            end
        end
        
        function val = tryGetPropertyValue(obj, name, zeroIsOK)
            %TRYGETPROPERTYVALUE returns the property value if it exists
            %   and is properly initialized (not zero or NaN). If ZEROISOK
            %   is set to TRUE, the property may have the value 0.
            %   Use it for non-required properties, which are required in
            %   certain functions (e.g. DENSITY for computeMassMatrix).
            if nargin == 2
                zeroIsOK = false;
            end
            
            if (~ obj.getProperties.hasValue(name)) ...
                    || isnan(obj.getPropertyValue(name))
                msg = [class(obj), ': ', name, ' in element ', ...
                    num2str(obj.getId()), ' is not initialized.'];
                e = MException('MATLAB:bm_mfem:missingProperty',msg);
                throw(e);
            elseif ~zeroIsOK && (obj.getPropertyValue(name) == 0)
                msg = [class(obj), ': ', name, ' in element ', ...
                    num2str(obj.getId()), ' is set to zero.'];
                e = MException('MATLAB:bm_mfem:missingProperty',msg);
                throw(e);
            else
                val = obj.getPropertyValue(name);
            end
            
        end
        
        function overwriteNode(element, oldNode, newNode)
            if ~isa(newNode,'Node')
                error('invalid node')
            end
            
            for itNode = 1:length(element.nodeArray)
                if (element.nodeArray(itNode) == oldNode)
                    
                    element.nodeArray(itNode) = newNode;
                    element.addDofsToSingleNode(newNode);
                    element.update;
                end
            end
        end
        
        function postprocess(obj)
        end
        
    end
    
    methods (Static, Sealed, Access = protected)
        function obj = getDefaultScalarElement
            obj = DefaultScalarElement;
        end
    end
    
end

