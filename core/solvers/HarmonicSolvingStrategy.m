classdef HarmonicSolvingStrategy < Solver
    %HARMONICSOLVINGSTRATEGY solves in case of prescribed dofs for
    % unkown dofs and forces;
    % If static solution is required set input parameter omega to zero;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  One can only put loads on unprescribed DOFS!!!!  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties (Access = private)
        freeDofs
        fixedDofs
        prescribedDofs
    end
    
    methods
        
        %constructor
        function obj = HarmonicSolvingStrategy(femModel, varargin)
            p = Solver.solverInputParser();
            p.parse(varargin{:});
            
            super_args = {femModel, p.Results.assembler, p.Results.rebuildMatrix};
            obj@Solver(super_args{:});
        end
        
        function solve(obj, omega)
            if ~obj.isInitialized; obj.initialize(); end
            
            if ~isempty(obj.prescribedDofs)
                if obj.rebuildMatrix; obj.buildFullSystem(); end
                x = obj.solvePrescribed(omega);
            else
                if obj.rebuildMatrix; obj.buildSystem(); end
                x = obj.solveHarmonic(omega);
            end
            
            SimpleAssembler.appendValuesToDofs(obj.femModel, x);
            obj.femModel.getProcessInfo.appendValue('FREQUENCY', omega);
        end
        
        function initialize(obj)
            obj.femModel.initialize();
            [obj.freeDofs, obj.fixedDofs, obj.prescribedDofs] = obj.femModel.getDofConstraints();
            
            if ~isempty(obj.prescribedDofs)
                %the prescribed solver needs the full matrices
                obj.buildFullSystem();
            else                
                obj.buildSystem();
            end
            
            obj.isInitialized = true;
            
        end
        
        function nodalForces = getNodalForces(obj)
            [ ~,nodalForces] = solve(obj);         
        end
        
    end
    
    methods (Access = private)
        function buildFullSystem(obj)
            %HARMONICSOLVINGSTRATEGY.BUILDFULLSYSTEM builds system and
            % does not include dirichlet conditions
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.dampingMatrix = obj.assembler.assembleGlobalDampingMatrix(obj.femModel);
            obj.forceVector = obj.assembler.applyExternalForces(obj.femModel);
        end
        
        function x = solveHarmonic(obj, omega)
            tmp = obj.stiffnessMatrix ...
                + 1i*omega*obj.dampingMatrix ...
                - omega^2*obj.massMatrix;
            x = tmp \ obj.forceVector;
        end
        
        function [xTotal, fTotal] = solvePrescribed(obj, omega)
            %HARMONICSOLVINGSTRATEGY.SOLVEPRESCRIBED solves for prescribed dofs
            %   xTotal: solved dof vector
            %   fTotal: solved force vector
            freeDofIds = obj.freeDofs.getId();
            prescribedDofIds = obj.prescribedDofs.getId();
            known_free = intersect(prescribedDofIds,freeDofIds);
            unknown_free = setdiff(freeDofIds,prescribedDofIds);
            
            [ k11, k12, k21, k22 ] = DisassembleMatrix(obj.stiffnessMatrix, unknown_free, known_free);
            [ m11, m12, m21, m22 ] = DisassembleMatrix(obj.massMatrix, unknown_free, known_free);
            [ d11, d12, d21, d22 ] = DisassembleMatrix(obj.dampingMatrix, unknown_free, known_free);
            [ fpre, ~ ] = DisassembleVector(obj.forceVector, unknown_free, known_free);
            
            upre = zeros(length(obj.prescribedDofs),1);
            for ii=1:length(obj.prescribedDofs)
                upre(ii)=obj.prescribedDofs(ii).getValue('end');
            end
            
            a11=k11+1i*omega*d11-omega^2*m11;
            a12=k12+1i*omega*d12-omega^2*m12;
            a21=k21+1i*omega*d21-omega^2*m21;
            a22=k22+1i*omega*d22-omega^2*m22;
            
            x = a11\(fpre - a12*upre);
            f = a21*x + a22*upre;
            
            dofs = obj.femModel.getDofArray;
            nDofs = length(dofs);
            
            xTotal=zeros(1,nDofs);
            xTotal(known_free)=upre;
            xTotal(unknown_free)=x;
            
            fTotal=zeros(1, nDofs);
            fTotal(known_free)=f;
            fTotal(unknown_free)=fpre;
            
            xTotal = applyVectorBoundaryConditions(xTotal, obj.fixedDofs.getId());
            
        end
    end
end
