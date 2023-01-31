classdef BlochInverse1DSolvingStrategy < Solver
    
    properties (Access = private)        
        Ksorted
        Msorted
        FixedDofs
        
        leftDofs
        rightDofs
        interiorDofs
        leftNodeIds
        rightNodeIds
    end
    
    methods
        function obj = BlochInverse1DSolvingStrategy(femModel, varargin)
            p = Solver.solverInputParser();
            p.parse(varargin{:});
            
            super_args = {femModel, p.Results.assembler, p.Results.rebuildMatrix};
            obj@Solver(super_args{:});
        end
        
        function [phases,frequencies] = solve(obj,numberOfPhases,numberOfBands)
            if ~ obj.isInitialized
                obj.initialize();
            end
            
            % calculate the transformed matrices
            [phases,mius] = obj.propConst(numberOfPhases);
            
            % initialize empty solution vector for frequencies
            % each row of represents one band
            frequencies=zeros(numberOfBands,length(phases));
            
            % iterate through all mius/phases
            for i=1:length(mius)
                miu=mius(i);
                
                % compute reduced matrices
                [Kred,Mred] = reducedStiffnesAndMass(obj,miu);
                
                % solve eigenvalue problem
                omega2 = eigs(Kred,Mred,numberOfBands,'sm');
                freq = sqrt(abs(omega2))/(2*pi);
                freq = sort(freq);
                
                % store the solution
                frequencies(:,i)=freq;
            end
        end
        
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end
            
            % find left and right Nodes
            nodeIdsLeft = obj.findLeftNodes();
            nodeIdsRight = obj.findRightNodes();
            
            % save it as property of the object
            obj.leftNodeIds = nodeIdsLeft;
            obj.rightNodeIds = nodeIdsRight;
            
            % find dof Ids of repective nodes
            [leftDofIds,rightDofIds,interiorDofIds] = getLeftRightInteriorDofIds(obj);
            
            % save it as property of the object
            obj.leftDofs = leftDofIds;
            obj.rightDofs = rightDofIds;
            obj.interiorDofs = interiorDofIds;
            
            % test the assignement of the nodes
            obj.testAssignmentOfNodes()
            
            % get FE-Matrices from femModel
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            
            % sort matrices with respect to left interior and right dof
            [obj.Ksorted,obj.Msorted] = obj.sortKandM();
             
        end %end initialize
        
        function testAssignmentOfNodes(obj)
            %TESTASSIGNMENTOFNODES tests if there is for every left node a
            %right node with the same y-coordinates.
            %Do left nodes have the same x-coordinate? Same for the right
            %nodes
 
            leftNodes = getNodes(obj.femModel, obj.leftNodeIds);
            rightNodes = getNodes(obj.femModel, obj.rightNodeIds);
            
            leftNodeX = getX(leftNodes);
            leftNodeY = getY(leftNodes);
            rightNodeX = getX(rightNodes);
            rightNodeY = getY(rightNodes);
            
            if length(obj.leftNodeIds) ~= length(obj.rightNodeIds)
                error('Same amount of left and right boundary nodes are requiered')
            end
            
            for i = 1:length(leftNodeY)
                if leftNodeY(i) ~= rightNodeY(i)
                    error('corresponding boundary nodes must have the same y-coordinates')
                end
            end
            
            for i = 1:length(leftNodeX)
                if leftNodeX(1) ~= leftNodeX(i)
                    error('All left boundary nodes must have the same x-coordinates')
                end
                if rightNodeX(1) ~= rightNodeX(i)
                    error('All right boundary nodes must have the same x-coordinates')
                end
            end
            
        end

        function [nodeIdsLeft] = findLeftNodes(obj)
            %FINDLEFTNODES gets the node-Ids of the left nodes and saves them
            %in the object properties
            nodeArray = obj.femModel.getAllNodes;
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            nodeXcoords = arrayfun(@(node) node.getX, nodeArray);
            
            sortedX = sort(nodeXcoords);
            minX = sortedX(1);
            l_xCoords = length(nodeXcoords);
            n=0;
            
            for i=1:l_xCoords
                if nodeXcoords(i) == minX
                    n = n+1;
                end
            end
            nodeIdsLeft = zeros(1,n);
            m = 0;
            for i=1:l_xCoords
                if nodeXcoords(i) == minX
                    m = m+1;
                    nodeIdsLeft(m) = nodeIds(i);
                end
            end
            %fprintf('Number of left boundary nodes is %s. \n', num2str(n))
            obj.leftNodeIds = nodeIdsLeft;   
        end

        function [nodeIdsRight] = findRightNodes(obj)
            %FINDRIGHTNODES gets the node-Ids of the right nodes and saves them
            %in the object properties
            nodeArray = obj.femModel.getAllNodes;
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            nodeXcoords = arrayfun(@(node) node.getX, nodeArray);
            sortedX = sort(nodeXcoords);
            
            maxX = sortedX(length(sortedX));
            l_xCoords = length(nodeXcoords);
            n=0;
            for i=1:l_xCoords
                if nodeXcoords(i) == maxX
                    n = n+1;
                end
            end
            m = 0;
            nodeIdsRight = zeros(1,n);
            for i=1:l_xCoords
                if nodeXcoords(i) == maxX
                    m = m+1;
                    nodeIdsRight(m) = nodeIds(i);
                end
            end
            obj.rightNodeIds = nodeIdsRight;
            %fprintf('Number of right boundary nodes is %s. \n', num2str(n))
        end
        
        function [leftDofIds,rightDofIds,interiorDofIds] = getLeftRightInteriorDofIds(obj)  %right and left nodes are requiered
            %GETLEFTRIGHTINTERIORDOFIDS gets the unfixed DOF IDs of the
            %left, right and interior nodes and saves them in the object
            %properties
            nodeArray = obj.femModel.getAllNodes;
            nodeIdsRight = obj.rightNodeIds;
            nodeIdsLeft = obj.leftNodeIds;
            
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            dofArray1 = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false)'; 
            
            indicesRightNodes=find(ismember(nodeIds,nodeIdsRight));
            indicesLeftNodes=find(ismember(nodeIds,nodeIdsLeft));
            
            dofArrayRight=dofArray1(indicesRightNodes); 
            dofArrayRight = [dofArrayRight{:}];
            
            dofArrayLeft=dofArray1(indicesLeftNodes);
            dofArrayLeft = [dofArrayLeft{:}];
            
            dofArrayInterior=dofArray1;
            dofArrayInterior([indicesRightNodes,indicesLeftNodes])=[];
            dofArrayInterior = [dofArrayInterior{:}];
            
            rightDofIds = arrayfun(@(dof) dof.getId, dofArrayRight);
            leftDofIds = arrayfun(@(dof) dof.getId, dofArrayLeft);
            interiorDofIds = arrayfun(@(dof) dof.getId, dofArrayInterior);
            
            %%Eliminate FixedDofs
            [~, fixedDofs] = getDofConstraints(obj.femModel);
            fixedDofIds = getId(fixedDofs);
            
            rightDofIds(any(rightDofIds==fixedDofIds'))=[];
            
            leftDofIds(any(leftDofIds==fixedDofIds'))=[];
            
            interiorDofIds(any(interiorDofIds==fixedDofIds'))=[];
%             indicesFixedInteriorIds = find(ismember(interiorDofIds,fixedDofIds));
%             interiorDofIds(indicesFixedInteriorIds)=[]; %same result, but
%             with a warning
            
            % save dofIds in object properties
            obj.leftDofs=leftDofIds;
            obj.rightDofs=rightDofIds;
            obj.interiorDofs=interiorDofIds;
            
        end
                
        function R = transformationMatrix(obj,miu)   
            %TRANSFORMATIONMATRIX creats transformationmatrix R, for a
            %specific miu
            numL = length(obj.leftDofs);
            numR = length(obj.rightDofs);
            numI = length(obj.interiorDofs);
            
            
            R = [eye(numL)       , zeros(numL,numI) ; ...
                zeros(numI,numL), eye(numI)        ; ...
                miu*eye(numR)   , zeros(numR,numI)];
        end
        
        function [Ksorted,Msorted] = sortKandM(obj)
            %SORTKANDM sorts stiffness and massmatrix, so that the entries
            %belonging to the respective DOFS of the left Nodes are on the
            %left-hand side, the right-DOF-entries are on the right-hand
            %side and the interior-DOF are in the middle
            vecdofsAll = [obj.leftDofs,obj.interiorDofs,obj.rightDofs];
            Ksorted = obj.stiffnessMatrix(vecdofsAll,vecdofsAll);
            Msorted = obj.massMatrix(vecdofsAll,vecdofsAll);
        end      
        
        function [Kred,Mred] = reducedStiffnesAndMass(obj,miu)
            %REDUCEDSTIFFNESSANDMASS gets the final stiffness and mass
            %matrices, which include the Bloch-Theorem
            R = transformationMatrix(obj,miu);
            
            Mred = R'*obj.Msorted*R;
            Kred = R'*obj.Ksorted*R;
        end

    end %end methods
    
    methods(Static)
        function [phase,miu] = propConst(numberOfPhases)
            %PROPCONST generates the phases (needed for the dispersion
            %curves) and the miu (needed for the transformationMatrix R); 
            phase = linspace(10e-3,pi,numberOfPhases); %10e-3 can be lowered up to 0, but then a warning appears 
            miu = exp(1i*phase); 
        end
        
    end %end static methods
    
    
    
end %end classdef


