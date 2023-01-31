classdef BlochDirect1DSolvingStrategy < Solver
    
    properties (Access = private)        
        Ksorted
        Msorted
        Dsorted
               
        leftDofs
        rightDofs
        interiorDofs
        leftNodes
        rightNodes
    end
    
    methods
        % Contructor
        function obj = BlochDirect1DSolvingStrategy(femModel, varargin)
            p = Solver.solverInputParser();
            p.parse(varargin{:});
            
            super_args = {femModel, p.Results.assembler, p.Results.rebuildMatrix};
            obj@Solver(super_args{:});
        end % end constructor
        
        % Solver function
        function [EV3]=solve(obj,frequencies)
        %function [EV1,EV2,EV3]=solve(obj,NumberOfFrequencies,maxFreq)   
            % initialize
            if ~ obj.isInitialized
                obj.initialize();
            end
            
            % transform frequency from Hz to rad/s
            omega = frequencies.*2*pi;
            tic            
            for k = 1:length(frequencies)
                omega_k = omega(k);
                % compute condensed dynamic matrix
                [Dcondensed] = condenseD(obj,omega_k);

                 % solve eigenvalue problem
%                     % option 1: Transfer Matrix T
%                     [T] = Tmatrix(obj,Dcondensed);
%                     EV = eig(full(T));
%                     EV1(:,k) = EV;
%                     % option 2: LR method
%                     [L,R] = LRmatrix(obj,Dcondensed);
%                     EV = eig(full(L),full(R));
%                     EV2(:,k) = EV;
                    % option 3: NL method
                    [N,L] = NLmatrix(obj,Dcondensed);
                    EV = eig(full(N),full(L));
                    EV3(:,k) = EV;
                    
                k = k+1;        
            end
                toc
        end % end solve
        
        % Methods
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end % end if
            
            % find left and right Nodes
            nodeIdsLeft = obj.findLeftNodes();
            nodeIdsRight = obj.findRightNodes();
            
            % save it as property of the object
            obj.leftNodes = nodeIdsLeft;
            obj.rightNodes = nodeIdsRight;
            
            % find dof Ids of repective nodes
            [leftDofIds,rightDofIds,interiorDofIds] = getLeftRightInteriorDofIds(obj);
            
            % save it as property of the object
            obj.leftDofs = leftDofIds;
            obj.rightDofs = rightDofIds;
            obj.interiorDofs = interiorDofIds;
            
            % test the assignement of the nodes
            %obj.testAssignmentOfNodes()
            
            % get FE-Matrices from femModel
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            obj.dampingMatrix = obj.assembler.assembleGlobalDampingMatrix(obj.femModel);

            % sort matrices with respect to left interior and right dof
            [obj.Ksorted,obj.Msorted,obj.Dsorted] = obj.sortKandM();
             
        end %end initialize
        
        function testAssignmentOfNodes(obj)
            leftNodes = getNodes(obj.femModel, obj.leftNodes);
            rightNodes = getNodes(obj.femModel, obj.rightNodes);
            
            leftNodeX = getX(leftNodes);
            leftNodeY = getY(leftNodes);
            rightNodeX = getX(rightNodes);
            rightNodeY = getY(rightNodes);
                        
            disp('left Nodes: [id,x,y]')
            X=[obj.leftNodes.' leftNodeX.' leftNodeY.'];
            disp(X)
            
            disp('right Nodes: [id,x,y]')
            Y=[obj.rightNodes.' rightNodeX.' rightNodeY.'];
            disp(Y) 
            
            if length(obj.leftNodes) ~= length(obj.rightNodes)
                error('Same amount of left and right boundary nodes are requiered')
            end % end if
            
            for i = 1:length(leftNodeY)
                if leftNodeY(i) ~= rightNodeY(i)
                    error('corresponding boundary nodes must have the same y-coordinates')
                end % end if
            end % end for
            
            for i = 1:length(leftNodeX)
                if leftNodeX(1) ~= leftNodeX(i)
                    error('All left boundary nodes must have the same x-coordinates')
                end % end if
                if rightNodeX(1) ~= rightNodeX(i)
                    error('All right boundary nodes must have the same x-coordinates')
                end % end if
            end % end for
            
        end % end testAssignmentOfNodes
        
        function [nodeIdsLeft] = findLeftNodes(obj)
            
            nodeArray = obj.femModel.getAllNodes;
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            nodeXcoords = arrayfun(@(node) node.getX, nodeArray);
            
            sortedX = sort(nodeXcoords);
            minX = sortedX(1);
            n=0;
            for i=1:length(nodeXcoords)
                if nodeXcoords(i) == minX
                    n = n+1;
                    nodeIdsLeft(n) = nodeIds(i);
                end % end if
            end % end for
            %fprintf('Number of left boundary nodes is %s. \n', num2str(n))
            obj.leftNodes = nodeIdsLeft;    %anstatt  obj.leftNodes(n) = nodeIds(i);
        end % end findLeftNodes
        
        function [nodeIdsRight] = findRightNodes(obj)
            
            nodeArray = obj.femModel.getAllNodes;
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            nodeXcoords = arrayfun(@(node) node.getX, nodeArray);
            sortedX = sort(nodeXcoords);
            
            maxX = sortedX(length(sortedX));
            n=0;
            for i=1:length(nodeXcoords)
                if nodeXcoords(i) == maxX
                    n = n+1;
                    nodeIdsRight(n) = nodeIds(i);
                end % end if
            end % end for
            obj.rightNodes = nodeIdsRight;
            %fprintf('Number of right boundary nodes is %s. \n', num2str(n))
        end % end findRightNodes
        
        function [leftDofIds,rightDofIds,interiorDofIds] = getLeftRightInteriorDofIds(obj)  %right and left nodes are requiered
            nodeArray = obj.femModel.getAllNodes;
            nodeIdsRight = obj.rightNodes;
            nodeIdsLeft = obj.leftNodes;
            
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
            
            % Eliminate FixedDofs
            
            [~, fixedDofs] = getDofConstraints(obj.femModel);
            fixedDofIds = getId(fixedDofs);
            
            rightDofIds(any(rightDofIds==fixedDofIds'))=[];
            
            leftDofIds(any(leftDofIds==fixedDofIds'))=[];
            
            interiorDofIds(any(interiorDofIds==fixedDofIds'))=[];
            
            % save dofIds in object properties
            obj.leftDofs=leftDofIds;
            obj.rightDofs=rightDofIds;
            obj.interiorDofs=interiorDofIds;
            
        end % end getLeftRightInteriorDofIds
        
        function [Ksorted,Msorted,Dsorted] = sortKandM(obj)
            vecdofsAll = [obj.leftDofs,obj.rightDofs,obj.interiorDofs];
            Ksorted = obj.stiffnessMatrix(vecdofsAll,vecdofsAll);
            Msorted = obj.massMatrix(vecdofsAll,vecdofsAll);
            Dsorted = obj.dampingMatrix(vecdofsAll,vecdofsAll);
        end % end sortKandM
        
        function [dynamicMatrix] = findD(obj,omega)
            dynamicMatrix = obj.Ksorted+1i*omega*obj.Dsorted-(omega^2).*obj.Msorted;
        end % end findD
        
        function [Dcondensed] = condenseD(obj,omega)
            
            dynamicMatrix = findD(obj,omega);
            
            nL = length(obj.leftDofs);
            nR = length(obj.rightDofs);
            nI = length(obj.interiorDofs);
            
            D_LL = dynamicMatrix(1:nL, 1:nL);
            D_RR = dynamicMatrix(nL+1:nL+nR, nL+1:nL+nR);
            D_II = dynamicMatrix(nL+nR+1:nL+nR+nI, nL+nR+1:nL+nR+nI);
            D_LR = dynamicMatrix(1:nL, nL+1:nL+nR);
            D_RL = D_LR.';
            D_LI = dynamicMatrix(1:nL, nL+nR+1:nL+nR+nI);
            D_IL = D_LI.';
            D_RI = dynamicMatrix(nL+1:nL+nR, nL+nR+1:nL+nR+nI);
            D_IR = D_RI.';
            
            D_LL_new = D_LL-D_LI*inv(D_II)*D_IL;
            D_RR_new = D_RR-D_RI*inv(D_II)*D_IR;
            D_LR_new = D_LR-D_LI*inv(D_II)*D_IR;
            %D_RL_new = D_RL-D_RI*inv(D_II)*D_IL;
            
            Dcondensed = [D_LL_new, D_LR_new; D_LR_new.', D_RR_new]; 
            
        end % end condenseD
        
        function [T] = Tmatrix(obj,Dcondensed)
            
            nL = length(obj.leftDofs);
            D_LL = Dcondensed(1:nL, 1:nL);
            D_RR = Dcondensed(nL+1:2*nL, nL+1:2*nL);
            D_LR = Dcondensed(1:nL, nL+1:2*nL);
            D_RL = D_LR.';
            
            T = [-inv(D_LR)*D_LL, inv(D_LR); ...
                -D_RL+D_RR*inv(D_LR)*D_LL, -D_RR*inv(D_LR)];
            
        end % end Tmatrix
        
        function [L,R] = LRmatrix(obj,Dcondensed)
               
            nL = length(obj.leftDofs);
            D_LL = Dcondensed(1:nL, 1:nL);
            D_RR = Dcondensed(nL+1:2*nL, nL+1:2*nL);
            D_LR = Dcondensed(1:nL, nL+1:2*nL);
            D_RL = D_LR.';
            
            L = [zeros(nL), D_RL; -D_RL, -D_LL-D_RR];
            R = [D_RL, zeros(nL); zeros(nL), D_LR];
            
        end % end LRmatrix
        
        function [N,L] = NLmatrix(obj,Dcondensed)
               
            nL = length(obj.leftDofs);
            D_LL = Dcondensed(1:nL, 1:nL);
            D_RR = Dcondensed(nL+1:2*nL, nL+1:2*nL);
            D_LR = Dcondensed(1:nL, nL+1:2*nL);
            D_RL = D_LR.';
            
            N = [zeros(nL), -eye(nL); -D_RL, D_RR];
            L = [eye(nL), zeros(nL); D_LL, -D_LR];
                        
        end % end NLmatrix
        
    end % end methods
    
end % classdef
