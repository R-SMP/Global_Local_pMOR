classdef Solver < handle
    %SOLVER Base class for solver implementations
    %   Detailed explanation goes here
    
    properties (Access = protected)
        femModel
        assembler
        
        stiffnessMatrix
        dampingMatrix
        massMatrix
        forceVector
        
        isInitialized
        rebuildMatrix
    end
    
    methods
        function obj = Solver(femModel, assembler, rebuildMatrixFlag)
            if nargin ~= 3
                msg = 'Wrong number of input arguments';
                e = MException('MATLAB:bm_mfem:SolverError',msg);
                throw(e);
            end
            
            obj.femModel = femModel;
            obj.assembler = assembler;
            obj.rebuildMatrix = rebuildMatrixFlag;
            
            obj.isInitialized = false;
        end
        
        function initialize(obj)
            obj.femModel.initialize();
            obj.buildSystem();
            obj.isInitialized = true;
        end
        
        function buildSystem(obj)
            [~, obj.stiffnessMatrix] = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            [~, obj.dampingMatrix] = obj.assembler.assembleGlobalDampingMatrix(obj.femModel);
            [~, obj.massMatrix] = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            [~, obj.forceVector] = obj.assembler.applyExternalForces(obj.femModel);
        end
    end
    
    methods (Sealed)
        function setIsInitialized(obj, flag)
            obj.isInitialized = flag;
        end
        
        function flag = getIsInitialized(obj)
            flag = obj.isInitialized;
        end
        
        function setRebuildMatrix(obj, flag)
            obj.rebuildMatrix = flag;
        end
        
        function flag = getRebuildMatrix(obj)
            flag = obj.rebuildMatrix;
        end
    end
    
    methods (Abstract)
        solve(obj)
    end
    
    methods (Static, Sealed)
        function p = solverInputParser()
            p = inputParser();
            p.addParameter('rebuildMatrix', false, @islogical);
            p.addParameter('assembler', SimpleAssembler, @(x) isa(x,'Assembler'));
        end
    end
    
end

