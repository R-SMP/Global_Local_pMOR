classdef Node < handle & matlab.mixin.Copyable
    %NODE The node class
    %   Parameters:
    %       id: unique identifier
    %       x, y(, z): coordinates
    
    properties (Access = private)
        id
        x
        y
        z
%         dofArray
        dofMap
        nProperties
        dimension
    end
    
    methods
        % constructor
        function node = Node(id, x, y, z)
            switch nargin
                case 0
                    % the empty constructor is needed in order to
                    % preallocate empty arrays of Nodes
                    node.id = -1;
                case 3
                    node.id = id;
                    node.x = x;
                    node.y = y;
                    node.dimension = 2;
                case 4
                    node.id = id;
                    node.x = x;
                    node.y = y;
                    node.z = z;
                    node.dimension = 3;
                otherwise
                    error('Wrong number of arguments')
            end
            node.dofMap = containers.Map;
            node.nProperties = PropertyContainer();
        end
        
        
        % getter functions
        function id = getId(node)
            %GETID Return the id of the node
            id = zeros;
            for ii = 1:length(node)
                id(ii) = node(ii).id;
            end
        end
        
        function coords = getCoords(nodes)
            %GETCOORDS Return all coordinates in the form x, y(, z)
            nnodes = length(nodes);
            coords = zeros(nnodes, nodes(1).dimension);
            for ii = 1:nnodes
                coords(ii,:) = [nodes(ii).x nodes(ii).y nodes(ii).z];
            end
        end
        
        function dim = getDimension(nodes)
            dim = zeros(1,length(nodes));
            for ii = 1:length(nodes)
                dim(ii) = nodes(ii).dimension;
            end
        end
        
        function x = getX(nodes)
            x = zeros;
            for ii = 1:length(nodes)
                x(ii) = nodes(ii).x;
            end
        end
        
        function y = getY(nodes)
            y = zeros;
            for ii = 1:length(nodes)
                y(ii) = nodes(ii).y;
            end
        end
        
        function z = getZ(nodes)
            z = zeros;
            for ii = 1:length(nodes)
                z(ii) = nodes(ii).z;
            end
        end
        
        function dofs = getDofArray(node)
            dofs = node.dofMap.values;
            dofs = [dofs{:}];
        end
        
        function dofs = getDofMap(node)
            dofs = node.dofMap;
        end

        function dof = getDof(nodes, dofName)
        %GETDOF get the dof with the specified name
        % parameters: node, dofName
            nNodes = length(nodes);
            dof = Dof.empty;
            for ii = 1:nNodes
                dof(ii) = nodes(ii).dofMap(dofName);
            end
        end        
        
        % member functions
        function addDof(nodes, dofNames)
            for itNode = 1:length(nodes)
                node = nodes(itNode);
                availableDofNames = node.dofMap.keys;
                for itDof = 1:length(dofNames)
                    name = dofNames{itDof};
                   if ~ any(ismember(availableDofNames,name))
                      newDof = Dof(node, 0.0, name);
                      node.dofMap(name) = newDof;
                   end
                end
            end
        end
        
        function removeDof(obj, dofName)
            for o = obj
                o.dofMap.remove(dofName);
            end
        end
        
        function fixDof(nodes, dofName)
            for ii = 1:length(nodes)
                dof = nodes(ii).dofMap(dofName);
                dof.fix();
            end
        end
        
        function fixAllDofs(nodes)
           for ii = 1:length(nodes)
              nodes(ii).getDofArray().fix(); 
           end
        end
        
        function unfixDof(nodes, dofName)
            for ii = 1:length(nodes)
                dof = nodes(ii).dofMap(dofName);
                dof.unfix();
            end
        end
        
        function prescribeDof(nodes, dofName, values)
            for ii = 1:length(nodes)
                dof = nodes(ii).dofMap(dofName);
                dof.prescribe(values(ii));
            end
        end
        
        function unprescribeDof(nodes, dofName)
            for ii = 1:length(nodes)
                dof = nodes(ii).dofMap(dofName);
                dof.unprescribe();
            end
        end
        
        function setDofValue(nodes, dofName, load, step)  %not load initial displacement
            %SETDOFVALUE set the value of a specific dof
            % parameters: dof, load
            if nargin == 3
                for ii = 1:length(nodes)
                    dof = nodes(ii).dofMap(dofName);
                    dof.setValue(load);
                end
            elseif nargin == 4
                for ii = 1:length(nodes)
                    dof = nodes(ii).dofMap(dofName);
                    dof.setValue(load,step);
                end
            else
                error('no')
            end
        end
        
        
        
        function setDofLoad(nodes, dofName, load)
            %SETDOFLOAD set the load of a specific dof
            % parameters: dof, load
            for ii = 1:length(nodes)
                dof = nodes(ii).dofMap(dofName);
                dof.setLoad(load);
            end
        end
        
        
        
        function val = getDofValue(obj, dofName, varargin)
            numvarargs = length(varargin);
            optargs = { 'end' };
            optargs(1:numvarargs) = varargin;
            step = optargs{:};
            
            nNodes = length(obj);
            
            if strcmp(step,'all')
                dof = obj(1).dofMap(dofName);
                tmp = dof.getValue(step);
                val = zeros(nNodes, length(tmp));
                val(1,:) = tmp;
                for ii = 2:nNodes
                    dof = obj(ii).dofMap(dofName);
                    val(ii,:) = dof.getValue(step);
                end
            else
                val = zeros(nNodes,1);
                for ii = 1:nNodes
                    if obj(ii).dofMap.isKey(dofName)
                        dof = obj(ii).dofMap(dofName);
                        val(ii) = dof.getValue(step);
                    else
                        val(ii) = 0;
                    end
                end
            end
        end
        
        function prop = getProperties(obj)
            prop = obj.nProperties;
        end
        
        function value = getPropertyValue(obj, valueName)
            for ii=1:length(obj)
                value(ii,:) = obj(ii).nProperties.getValue(valueName)';
            end
        end
        
        function setProperties(obj, props)
            for ii = 1:length(obj)
                names = props.getValueNames();
                for jj = 1:length(names)
                    obj(ii).nProperties.setValue(names{jj}, props.getValue(names{jj}));
                end
            end
        end
        
        function setPropertyValue(obj, valueName, value)
            for ii = 1:length(obj)
                obj(ii).nProperties.setValue(valueName, value);
            end
        end
        
        function addProperty(obj, valueName, value)
            for ii = 1:length(obj)
                obj(ii).nProperties.addValue(valueName, value);
            end
        end
        
        function addNewValue(nodes, valueNames)
            warning('deprecated function Node.addNewValue()')
            for itName = 1:length(valueNames)
                for itNode = 1:length(nodes)
                    nodes(itNode).nProperties.addValue(char(valueNames(itName)));
                end
            end
        end
        
        function setStepValue(nodes, valueName, value, step)
           for ii = 1:length(nodes)
              nodes(ii).nProperties.setStepValue(valueName, value, step); 
           end
        end
        
        function appendStepValue(nodes, valueName, value)
           for ii = 1:length(nodes)
              nodes(ii).nProperties.appendStepValue(valueName, value); 
           end
        end
        
        function val = getValue(nodes, valueName, step)
            warning('deprecated function Node.getValue')
            nNodes = length(nodes);
            type = checkPropertyName(valueName);
            if nargin == 2
                if (strcmp(type,'variable3d'))
                    error('Error: please specify a step to retrieve 3d values')
                end
%                 val = zeros(1,nNodes);
                for ii = 1:nNodes
                    val(ii,:) = nodes(ii).nProperties.getValue(valueName);
                end
            elseif nargin == 3
%                 if (strcmp(type,'variable3d'))
%                     val = zeros
                 for ii = 1:nNodes
                    val(ii,:) = nodes(ii).nProperties.getValue(valueName, step);
                end
            end
        end        
        
    end
    
    methods (Access = protected)
        
        function cp = copyElement(obj)
           % copy constructor
           %cp = copyElement@matlab.mixin.Copyable(obj);
            coords = obj.getCoords;
            if (length(coords) == 2)
                cp = Node(obj.getId, coords(1), coords(2));
            else
                cp = Node(obj.getId, coords(1), coords(2), coords(3));
            end
        end
    end
    
end

