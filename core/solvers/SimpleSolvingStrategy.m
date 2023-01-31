classdef SimpleSolvingStrategy < Solver
    %SIMPLESOLVINGSTRATEGY Solves for the static displacement
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        %constructor
        function obj = SimpleSolvingStrategy(femModel, varargin)
            p = Solver.solverInputParser();
            p.parse(varargin{:});
            
            super_args = {femModel, p.Results.assembler, p.Results.rebuildMatrix};
            obj@Solver(super_args{:});
        end
        
        function x = solve(obj)
            if ~ obj.isInitialized
                obj.initialize();
            end
            
            [~, Kred] = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            [~, fred] = obj.assembler.applyExternalForces(obj.femModel);
            
            x = Kred \ fred;
            
            obj.assembler.assignResultsToDofs(obj.femModel, x);
        end
        
        %SIMPLESOLVINGSTRATEGY.INITIALIZE overwrite base class, because the
        % stiffness matrix is only needed once
        function initialize(obj)
            obj.femModel.initialize;
            obj.isInitialized = true;
        end
        
        function nodalForces = getNodalForces(obj, step)
            nodalForces = SimpleAssembler.applyExternalForces(obj.femModel);
            [~, fixedDofs] = obj.femModel.getDofConstraints();
            
            K = SimpleAssembler.assembleGlobalStiffnessMatrix(obj.femModel);
            nodalForces(fixedDofs.getId) = K(fixedDofs.getId, :) ...
                * obj.femModel.getDofArray.getValue(step);
            
        end
        
    end
end
