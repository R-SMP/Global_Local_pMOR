classdef PorousBlochDirect1DSolvingStrategy < Solver
    
    properties (Access = private)
        dynamicMatrix
        
        Ksorted
        Msorted
               
        leftDofs
        rightDofs
        interiorDofs
        leftNodes
        rightNodes
        sortedDofs
        
        EigValue
        EigVector
        frequency
        k_I
        k_R
        cleanEigValue
    end
    
    methods
        % Contructor
        function obj = PorousBlochDirect1DSolvingStrategy(femModel, varargin)
            p = Solver.solverInputParser();
            p.parse(varargin{:});
            
            super_args = {femModel, p.Results.assembler, p.Results.rebuildMatrix};
            obj@Solver(super_args{:});
        end % end constructor
        
        % Solver function
        function [cleanEV, EVec, k_I, k_R] = solve(obj,omega_k,Lx)

            % initialize
            if ~ obj.isInitialized
                obj.initialize();
            end
%          n_ii = log10(abs(obj.Ksorted)./abs(omega_k^2*obj.Msorted));
%          n_ii(isinf(n_ii)) = 0;
%          n_ii(isnan(n_ii)) = 0;
%          n_max = max(max(n_ii));
            % compute condensed dynamic matrix
            [Dcondensed] = condenseD(obj,omega_k);

             % solve eigenvalue problem
%                     % option 1: Transfer Matrix T
%                     [T] = Tmatrix(obj,Dcondensed);
%                     EV = eig(full(T));
%                     
%                     % option 2: LR method
%                     [L,R] = LRmatrix(obj,Dcondensed);
%                     [EV,ev] = eig(full(L),full(R));
%                     
                % option 3: NL method
                [N,L] = NLmatrix(obj,Dcondensed);
                % Calculate Eigenvectors and Eigenvalues
                [EV,ev] = eig(full(N),full(L));
                
                % Write Eigenvalues in a Vector
                [obj.EigValue] = ev*ones(2*length(obj.leftDofs),1);
                % Decondense the inner Nodes and sort w.r.t Dofs
                [EV2] = decondenseEV(obj,EV);
                [obj.EigVector] = sortEigVectorsByDofs(obj,EV2);
                EVec = obj.EigVector;
                [k_I,k_R,cleanEV] = calculateWavenumbers(obj,obj.EigValue,Lx);
                
                
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
            obj.testAssignmentOfNodes()
            
            % get FE-Matrices from femModel
            [,obj.massMatrix] = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            [,obj.stiffnessMatrix] = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);

            % sort matrices with respect to left interior and right dof
            [obj.Ksorted,obj.Msorted] = obj.sortKandM();
             
        end %end initialize
        
        function testAssignmentOfNodes(obj)
            leftNodes = getNodes(obj.femModel, obj.leftNodes);
            rightNodes = getNodes(obj.femModel, obj.rightNodes);
            
            leftNodeX = getX(leftNodes);
            leftNodeY = getY(leftNodes);
            rightNodeX = getX(rightNodes);
            rightNodeY = getY(rightNodes);
                        
%             disp('left Nodes: [id,x,y]')
            X=[obj.leftNodes.' leftNodeX.' leftNodeY.'];
%             disp(X)
            
%             disp('right Nodes: [id,x,y]')
            Y=[obj.rightNodes.' rightNodeX.' rightNodeY.'];
%             disp(Y) 
            
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
%             fprintf('Number of left boundary nodes is %s. \n', num2str(n))
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
%             fprintf('Number of right boundary nodes is %s. \n', num2str(n))
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
        
        function [Ksorted,Msorted] = sortKandM(obj)
            vecdofsAll = [obj.leftDofs,obj.rightDofs,obj.interiorDofs];
            obj.sortedDofs = vecdofsAll;
            Ksorted = obj.stiffnessMatrix(vecdofsAll,vecdofsAll);
            Msorted = obj.massMatrix(vecdofsAll,vecdofsAll);
        end % end sortKandM
        
        function [dynamicMatrix] = findD(obj,omega)
            dynamicMatrix = obj.Ksorted-(omega^2).*obj.Msorted;
        end % end findD
        
        function [Dcondensed] = condenseD(obj,omega)
            
            dynamicMatrix = findD(obj,omega);
            obj.dynamicMatrix = dynamicMatrix;
            
            nL = length(obj.leftDofs);
            nR = length(obj.rightDofs);
            nI = length(obj.interiorDofs);
            
            D_LL = dynamicMatrix(1:nL, 1:nL);
            D_RR = dynamicMatrix(nL+1:nL+nR, nL+1:nL+nR);
            D_II = dynamicMatrix(nL+nR+1:nL+nR+nI, nL+nR+1:nL+nR+nI);
            D_LR = dynamicMatrix(1:nL, nL+1:nL+nR);
            D_RL = dynamicMatrix(nL+1:nL+nR, 1:nL);
            D_LI = dynamicMatrix(1:nL, nL+nR+1:nL+nR+nI);
            D_IL = dynamicMatrix(nL+nR+1:nL+nR+nI, 1:nL);
            D_RI = dynamicMatrix(nL+1:nL+nR, nL+nR+1:nL+nR+nI);
            D_IR = dynamicMatrix(nL+nR+1:nL+nR+nI, nL+1:nL+nR);
            
            D_LL_new = D_LL-D_LI*inv(D_II)*D_IL;
            D_RR_new = D_RR-D_RI*inv(D_II)*D_IR;
            D_LR_new = D_LR-D_LI*inv(D_II)*D_IR;
            D_RL_new = D_RL-D_RI*inv(D_II)*D_IL;
            
            Dcondensed = [D_LL_new, D_LR_new; D_RL_new, D_RR_new]; 
            
        end % end condenseD
        
        function [fullEV] = decondenseEV(obj,EV)
             
            nL = length(obj.leftDofs);
            nR = length(obj.rightDofs);
            nI = length(obj.interiorDofs);
            
            D_II = obj.dynamicMatrix(nL+nR+1:nL+nR+nI, nL+nR+1:nL+nR+nI);
            inv_D_II = inv(D_II);
            D_LI = obj.dynamicMatrix(1:nL, nL+nR+1:nL+nR+nI);
            D_IL = obj.dynamicMatrix(nL+nR+1:nL+nR+nI, 1:nL);
            D_RI = obj.dynamicMatrix(nL+1:nL+nR, nL+nR+1:nL+nR+nI);
            D_IR = obj.dynamicMatrix(nL+nR+1:nL+nR+nI, nL+1:nL+nR);
            
            U = [eye(nL),        zeros(nL,nL);
                 zeros(nL,nL),   -eye(nL);
                 -inv_D_II*D_IL, inv_D_II*D_IR];
            
            fullEV= U*EV;
            
        end % end decondenseEV
        
        function [sortedEV] = sortEigVectorsByDofs(obj,EV)
            % Sort EigVectors by Dofs
            sortedDofs = obj.sortedDofs';
            S = [sortedDofs, EV(:,:)];
            SortDofs = sortrows(S);
            sortedEV = SortDofs(:,2:size(S,2)); 
            
        end % end sortEigVectorsByDofs
        
        function [T] = Tmatrix(obj,Dcondensed)
            
            nL = length(obj.leftDofs);
            D_LL = Dcondensed(1:nL, 1:nL);
            D_RR = Dcondensed(nL+1:2*nL, nL+1:2*nL);
            D_LR = Dcondensed(1:nL, nL+1:2*nL);
            D_RL = Dcondensed(nL+1:2*nL, 1:nL);
            
            T = [-inv(D_LR)*D_LL, inv(D_LR); ...
                -D_RL+D_RR*inv(D_LR)*D_LL, -D_RR*inv(D_LR)];
            
        end % end Tmatrix
        
        function [L,R] = LRmatrix(obj,Dcondensed)
               
            nL = length(obj.leftDofs);
            D_LL = Dcondensed(1:nL, 1:nL);
            D_RR = Dcondensed(nL+1:2*nL, nL+1:2*nL);
            D_LR = Dcondensed(1:nL, nL+1:2*nL);
            D_RL = Dcondensed(nL+1:2*nL, 1:nL);
            
            L = [zeros(nL), D_RL; -D_RL, -D_LL-D_RR];
            R = [D_RL, zeros(nL); zeros(nL), D_LR];
            
        end % end LRmatrix
        
        function [N,L] = NLmatrix(obj,Dcondensed)
               
            nL = length(obj.leftDofs);
            D_LL = Dcondensed(1:nL, 1:nL);
            D_RR = Dcondensed(nL+1:2*nL, nL+1:2*nL);
            D_LR = Dcondensed(1:nL, nL+1:2*nL);
            D_RL = Dcondensed(nL+1:2*nL, 1:nL);
            
            N = [zeros(nL), -eye(nL); -D_RL, D_RR];
            L = [eye(nL), zeros(nL); D_LL, -D_LR];
                        
        end % end NLmatrix
        
        function postprocessingEV(obj,EigValue,EigVector)
            
            SimpleAssembler.appendValuesToDofs(obj.femModel, EigVector);
            obj.femModel.getProcessInfo.appendValue('FREQUENCY', EigValue);
            
        end % end postprocessingEV
        
        function [k_I,k_R,EV] = calculateWavenumbers(obj,EV,Lx)
            
            [EV] = cleanEigValues(obj,EV);
            
            lnEV = log(EV);
            k_I = real(lnEV)./Lx;
            k_R = -imag(lnEV);
            
            [k_I,k_R] = cleanWavenumbers(obj,k_I,k_R);
            
            obj.k_I = k_I;
            obj.k_R = k_R;
            
        end % end calculateWavenumbers
        
        function [EV] = cleanEigValues(obj,EV)
            nEV = size(EV,1);
            % Clean Noise 
            EV = EV;
            for s = 1:nEV
                if abs(EV(s))<1e-4
                   EV(s)=0;
                end
                if abs(EV(s))>1e3
                   EV(s)=0;
                end
            end
            obj.cleanEigValue = EV;
        end % cleanEigValues
        
        function [k_I,k_R] = cleanWavenumbers(obj,k_I,k_R)
            nEV = size(k_I,1);
            % Dismiss results with very high damping            
            for s = 1:nEV
                if k_I(s)<-1
                    k_I(s)=0;
                    k_R(s)=0;
                end
                if k_I(s)>1
                    k_I(s)=0;
                    k_R(s)=0;
                end
            end
            
        end % end cleanWavenumbers   
        
        
    end % end methods
    
end % classdef
