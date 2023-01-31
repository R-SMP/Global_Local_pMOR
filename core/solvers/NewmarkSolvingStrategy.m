classdef NewmarkSolvingStrategy < Solver
    %NEWMARKSOLVINGSTRATEGY A time integration method based on the Newmark
    %method
    %   Newmark's method as outlined in G�radin, Michel, and Daniel J. 
    %   Rixen. Mechanical vibrations: theory and application to structural
    %   dynamics. John Wiley & Sons, 2014.
    
    properties (Access = private)
        dt
        alpha
        beta
        gamma
        newmarkCoefficients
        
        lhs
    end
    
    methods
        
        %constructor
        function obj = NewmarkSolvingStrategy(femModel, dt, varargin)
            p = Solver.solverInputParser();
            p.addParameter('alpha', 0, @isnumeric);
            p.addParameter('beta', .25, @isnumeric);
            p.addParameter('gamma', .5, @isnumeric);
            p.parse(varargin{:});
            
            super_args = {femModel, p.Results.assembler, p.Results.rebuildMatrix};
            obj@Solver(super_args{:});
            
            obj.dt = dt;
            
            obj.alpha = p.Results.alpha;
            obj.beta = p.Results.beta;
            obj.gamma = p.Results.gamma;
        end
        
        function solve(solver)
            if ~solver.isInitialized
                solver.initialize;
            end
            step = solver.femModel.getProcessInfo.getValue('STEP');
            [~, fixedDofs] = solver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            [~, force] = solver.assembler.applyExternalForces(solver.femModel);
            dispOld = solver.assembler.assembleValuesVector(solver.femModel, step);
            dispOld = applyVectorBoundaryConditions(dispOld, fixedDofIds)';
            velOld = solver.assembler.assembleFirstDerivativesVector(solver.femModel, step);
%             velOld = solver.assembler.assembleFirstDerivativesVector2(solver.femModel, 'DISPLACEMENT_X', step) ...
%                 + solver.assembler.assembleFirstDerivativesVector2(solver.femModel, 'DISPLACEMENT_Y', step) ...
%                 + solver.assembler.assembleFirstDerivativesVector2(solver.femModel, 'DISPLACEMENT_Z', step);
            velOld = applyVectorBoundaryConditions(velOld, fixedDofIds)';
            accOld = solver.assembler.assembleSecondDerivativesVector(solver.femModel, step);
%             accOld = solver.assembler.assembleSecondDerivativesVector2(solver.femModel, 'DISPLACEMENT_X', step) ...
%                 + solver.assembler.assembleSecondDerivativesVector2(solver.femModel, 'DISPLACEMENT_Y', step) ...
%                 + solver.assembler.assembleSecondDerivativesVector2(solver.femModel, 'DISPLACEMENT_Z', step);
            accOld = applyVectorBoundaryConditions(accOld, fixedDofIds)';
            
            if solver.alpha == 0
                rhs = force ...
                    - solver.dampingMatrix * (velOld + ((1 - solver.gamma) .* solver.dt) * accOld) ...
                    - solver.stiffnessMatrix * (dispOld + solver.dt .* velOld + ((0.5 - solver.beta) * power(solver.dt,2)) .* accOld);
                
                accNew = solver.lhs \ rhs;
                velNew = velOld + ((1 - solver.gamma) * solver.dt) .* accOld + (solver.gamma * solver.dt) .* accNew;
                dispNew = dispOld + solver.dt .* velOld + (power(solver.dt, 2) * (0.5 - solver.beta)) .* accOld + (power(solver.dt, 2) * solver.beta) .* accNew;
                
                
            else
                c = solver.newmarkCoefficients;
                rhs = force ...
                    + solver.massMatrix * ((1 - solver.alpha) .* (dispOld .* c(1) + velOld .* c(3) + accOld .* c(4)) - solver.alpha .* accOld) ...
                    + solver.dampingMatrix * (dispOld .* c(2) + velOld .* c(5) + accOld .* c(6));
                
                dispNew = solver.lhs \ rhs;
                accNew = (dispNew - dispOld) .* c(1) - velOld .* c(3) - accOld .* c(4);
                velNew = velOld + (solver.dt * (1 - solver.gamma)) .* accOld + accNew .* (solver.dt * solver.gamma);
            end
            
            %write the values back to the dofs / nodes
            solver.assembler.appendValuesToDofs(solver.femModel, dispNew);
            solver.assembler.appendFirstDerivativeValuesToDofs(solver.femModel, velNew);
            solver.assembler.appendSecondDerivativeValuesToDofs(solver.femModel, accNew);
            
            solver.femModel.getProcessInfo.setValue('STEP', step+1);
            solver.femModel.getProcessInfo.appendValue('TIME', step*solver.dt);
            
        end
        
        function initialize(solver)
            solver.femModel.initialize;
            step = solver.femModel.getProcessInfo.getValue('STEP');
            
            % assemble and reduce matrices
            [~, fixedDofs] = solver.femModel.getDofConstraints();
            fixedDofIds = fixedDofs.getId();
            
            mass = solver.assembler.assembleGlobalMassMatrix(solver.femModel);
            solver.massMatrix = applyMatrixBoundaryConditions(mass, fixedDofIds);
            damping = solver.assembler.assembleGlobalDampingMatrix(solver.femModel);
            solver.dampingMatrix = applyMatrixBoundaryConditions(damping, fixedDofIds);
            stiffness = solver.assembler.assembleGlobalStiffnessMatrix(solver.femModel);
            solver.stiffnessMatrix = applyMatrixBoundaryConditions(stiffness, fixedDofIds);
            
            % initial acceleration
            [~, force0] = solver.assembler.applyExternalForces(solver.femModel);
            disp = solver.assembler.assembleValuesVector(solver.femModel, step);
            disp0 = applyVectorBoundaryConditions(disp, fixedDofIds)';
            vel = solver.assembler.assembleFirstDerivativesVector(solver.femModel, step);
            vel0 = applyVectorBoundaryConditions(vel, fixedDofIds)';
            
%             acc0 = (solver.massMatrix) \ (force0' - solver.stiffnessMatrix * disp0 - solver.dampingMatrix * vel0);
            acc0 = solver.massMatrix \ (force0 - solver.stiffnessMatrix * disp0 - solver.dampingMatrix * vel0);
            solver.assembler.setSecondDerivativeValuesToDofs(solver.femModel, acc0);
            solver.isInitialized = true;
            
            if solver.alpha == 0
                %set up left hand side of the les
                solver.lhs = solver.massMatrix ...
                    + (solver.gamma * solver.dt) .* solver.dampingMatrix ...
                    + (solver.beta * power(solver.dt,2)) .* solver.stiffnessMatrix;
            else
                %coefficients for bossak/newmark
                solver.beta = power((1.0 - solver.alpha),2) * solver.beta;
                solver.gamma = solver.gamma - solver.alpha;
                c = zeros(6,1);
                c(1) = 1 / (solver.beta * power(solver.dt,2));
                c(2) = solver.gamma / (solver.beta * solver.dt);
                c(3) = 1 / (solver.beta * solver.dt);
                c(4) = 1 / (2 * solver.beta) - 1;
                c(5) = solver.gamma / solver.beta - 1;
                c(6) = solver.dt * 0.5 * (solver.gamma / solver.beta - 2);
                solver.newmarkCoefficients = c;
                
                solver.lhs = solver.stiffnessMatrix ...
                    + ((1 - solver.alpha) * c(1)) .* solver.massMatrix ...
                    + c(2) .* solver.dampingMatrix;
            end
        end
    end
    
end

