classdef ApproximateClassicalElement2d4n < PorousElement2d4n
    %   ApproximateClassicalPorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Fluid-Displacement (see: PANNETON 1996/ RUMPLER 2013)
    %   Approximate Version with frequency independent matrices
    
    %   Coupling between a ClassicalPorousElement and a structure is
    %   considered here, but yet works only for 2-node interfaces between elements 
    %   belonging to different media (not 3 or 1 node interfaces)

    
    properties (Access = private)
    end
    
    methods
        
        % constructor
        function obj = ApproximateClassicalElement2d4n(id,nodeArray)
            requiredPropertyNames = cellstr(["DENSITY_SOLID", "LAMBDA_SOLID", ...
                "MUE_SOLID", "DAMPING_SOLID", "DENSITY_FLUID", ...
                "VISCOSITY_FLUID", "STANDARD_PRESSURE_FLUID", ...
                "HEAT_CAPACITY_FLUID", "PRANDTL_NUMBER_FLUID", "POROSITY", ...
                "TORTUOSITY", "FLOW_RESISTIVITY", "VISCOUS_LENGHT", ...
                "THERMAL_LENGTH", "NUMBER_GAUSS_POINT"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class contructor
            obj@PorousElement2d4n(super_args{:});
            obj.dofNames = cellstr(["DISPLACEMENT_X"; "DISPLACEMENT_Y"; ...
                "DISPLACEMENT_FLUID_X"; "DISPLACEMENT_FLUID_Y"]);
        end
        
        %Initialization
        function initialize(obj)
            obj.lengthX = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(2).getCoords);
            
            obj.lengthY = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(4).getCoords);
            
            checkConvexity(obj);
            
            % Default setting for parameters of fluid (AIR) and NUMBER_GAUSS_POINT
            obj.eProperties.setValue('DENSITY_FLUID',1.21);
            obj.eProperties.setValue('VISCOSITY_FLUID',1.84e-5);
            obj.eProperties.setValue('STANDARD_PRESSURE_FLUID',101e3);
            obj.eProperties.setValue('HEAT_CAPACITY_FLUID',1.4);
            obj.eProperties.setValue('PRANDTL_NUMBER_FLUID',0.71);
            obj.eProperties.setValue('NUMBER_GAUSS_POINT',2);
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            % Read Material Properties
            LAMBDA_SOLID = obj.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = obj.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = obj.getPropertyValue('DAMPING_SOLID');
            STANDARD_PRESSURE_FLUID = obj.getPropertyValue('STANDARD_PRESSURE_FLUID');
            POROSITY = obj.getPropertyValue('POROSITY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');

            % Determine Coefficients for EMat
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            
            D_1 = [1 0 0 ; 0 1 0 ; 0 0 0.5];
            D_2 = [1 1 0 ; 1 1 0 ; 0 0 0];
            
            D_S = 2*MUE*D_1 + (LAMBDA + ((1 - POROSITY)^2/POROSITY)*STANDARD_PRESSURE_FLUID)*D_2;
            D_SF = (1 - POROSITY)*STANDARD_PRESSURE_FLUID*D_2;
            D_F = POROSITY*STANDARD_PRESSURE_FLUID*D_2;
            
            stiffnessMatrix_S = sparse(8,8);
            stiffnessMatrix_SF = sparse(8,8);
            stiffnessMatrix_F = sparse(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, ~, Be, ~, J] = computeShapeFunction(obj, xi, eta);
                    stiffnessMatrix_S=stiffnessMatrix_S+(w(i)*w(j)*det(J)*transpose(Be)*(D_S*Be));
                    stiffnessMatrix_F=stiffnessMatrix_F+(w(i)*w(j)*det(J)*transpose(Be)*(D_F*Be));
                    stiffnessMatrix_SF=stiffnessMatrix_SF+(w(i)*w(j)*det(J)*transpose(Be)*(D_SF*Be));
                end
            end
            
            stiffnessMatrix = [stiffnessMatrix_S stiffnessMatrix_SF; stiffnessMatrix_SF stiffnessMatrix_F];    
            
            for i=1:4
                check(i) = obj.nodeArray(i).getPropertyValue('IS_COUPLING_NODE');
            end
            
            if sum(check)==2
                couplingIds = find(check);
                couplingNodes = obj.nodeArray(couplingIds);
                matrixIds = (couplingIds -1)*2 +1;
                normalVector = obj.getNormalVector(couplingNodes);
                
                if abs(normalVector) == [1 0]
                   stiffnessMatrix(:,matrixIds) = stiffnessMatrix(:,matrixIds) + stiffnessMatrix(:,matrixIds+8);
                   stiffnessMatrix(matrixIds,:) = stiffnessMatrix(matrixIds,:) + stiffnessMatrix(matrixIds+8,:);
                   
                   stiffnessMatrix(:,matrixIds+8) = 0;
                   stiffnessMatrix(matrixIds+8,:) = 0;
                   stiffnessMatrix(matrixIds(1)+8,matrixIds(1)+8) = 1;
                   stiffnessMatrix(matrixIds(2)+8,matrixIds(2)+8) = 1;
                   stiffnessMatrix(matrixIds(1),matrixIds(1)+8) = -1;
                   stiffnessMatrix(matrixIds(1)+8,matrixIds(1)) = -1;
                   stiffnessMatrix(matrixIds(2),matrixIds(2)+8) = -1;
                   stiffnessMatrix(matrixIds(2)+8,matrixIds(2)) = -1;
                          
                elseif abs(normalVector) == [0 1]
                   stiffnessMatrix(:,matrixIds+1) = stiffnessMatrix(:,matrixIds+1) + stiffnessMatrix(:,matrixIds+9);
                   stiffnessMatrix(matrixIds+1,:) = stiffnessMatrix(matrixIds+1,:) + stiffnessMatrix(matrixIds+9,:);
                   
                   stiffnessMatrix(:,matrixIds+9) = 0;
                   stiffnessMatrix(matrixIds+9,:) = 0;
                   stiffnessMatrix(matrixIds(1)+9,matrixIds(1)+9) = 1;
                   stiffnessMatrix(matrixIds(2)+9,matrixIds(2)+9) = 1;
                   stiffnessMatrix(matrixIds(1)+1,matrixIds(1)+9) = -1;
                   stiffnessMatrix(matrixIds(1)+9,matrixIds(1)+1) = -1;
                   stiffnessMatrix(matrixIds(2)+1,matrixIds(2)+9) = -1;
                   stiffnessMatrix(matrixIds(2)+9,matrixIds(2)+1) = -1;
                    
                else
                    temp1 = stiffnessMatrix(:,matrixIds+8)*normalVector(1);
                    temp2 = stiffnessMatrix(:,matrixIds+9)*normalVector(2);
                    
                    stiffnessMatrix(:,matrixIds) = stiffnessMatrix(:,matrixIds) + temp1;
                    stiffnessMatrix(:,matrixIds+1) = stiffnessMatrix(:,matrixIds+1) + temp2;
                    
                    stiffnessMatrix(matrixIds,:) = stiffnessMatrix(matrixIds,:) + temp1.';
                    stiffnessMatrix(matrixIds+1,:) = stiffnessMatrix(matrixIds+1,:) + temp2.';
                    
                    stiffnessMatrix(:,matrixIds+8) = stiffnessMatrix(:,matrixIds+8) - temp1;
                    stiffnessMatrix(:,matrixIds+9) = stiffnessMatrix(:,matrixIds+9) - temp2;
                    
                    stiffnessMatrix(matrixIds+8,:) = stiffnessMatrix(matrixIds+8,:) - temp1.';
                    stiffnessMatrix(matrixIds+9,:) = stiffnessMatrix(matrixIds+9,:) - temp2.';
                    
                end 
            elseif sum(check)==1 || sum(check)==3 || sum(check)==4
                e = MException('MATLAB:bm_mfem:cornerCoupling','Code can handle only 2-node coupling interfaces');
                throw(e);
            end         
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            % Read Material Properties
            DENSITY_SOLID = obj.getPropertyValue('DENSITY_SOLID'); 
            DENSITY_FLUID = obj.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = obj.getPropertyValue('VISCOSITY_FLUID');
            POROSITY = obj.getPropertyValue('POROSITY');
            TORTUOSITY = obj.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = obj.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = obj.getPropertyValue('VISCOUS_LENGHT');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine relevant densities
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            H = FLOW_RESISTIVITY^2*VISCOUS_LENGHT^2*POROSITY^2/(4*TORTUOSITY^2*VISCOSITY_FLUID*DENSITY_FLUID);
            alpha = POROSITY^2*FLOW_RESISTIVITY/(2*H);
            
            ElementmassMatrix=sparse(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(obj, xi, eta);
                    ElementmassMatrix = ElementmassMatrix + (w(i) * w(j) * transpose(N_mat) * N_mat * det(J));
                end
            end
            
            M = [((1-POROSITY)*DENSITY_SOLID + APPARENT_MASS_DENSITY) * ElementmassMatrix, - APPARENT_MASS_DENSITY * ElementmassMatrix; ...
                - APPARENT_MASS_DENSITY * ElementmassMatrix, (POROSITY*DENSITY_FLUID + APPARENT_MASS_DENSITY) * ElementmassMatrix];
        
            ElementDampingMatrix=sparse(8,8);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(obj, xi, eta);
                    ElementDampingMatrix = ElementDampingMatrix + (w(i) * w(j) * transpose(N_mat) * N_mat * det(J));
                end
            end
            
            C = [ElementDampingMatrix, - ElementDampingMatrix; ...
                - ElementDampingMatrix, ElementDampingMatrix];
            
            massMatrix = M - i*alpha*C;
            massMatrix = obj.PorousStructureCoupling(massMatrix);
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            % Read Material Properties
            FLOW_RESISTIVITY = obj.getPropertyValue('FLOW_RESISTIVITY');
            POROSITY = obj.getPropertyValue('POROSITY');
            DENSITY_FLUID = obj.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = obj.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = obj.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = obj.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = obj.getPropertyValue('PRANDTL_NUMBER_FLUID');
            THERMAL_LENGTH = obj.getPropertyValue('THERMAL_LENGTH');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            beta = POROSITY^2 * FLOW_RESISTIVITY;
            H = 16 * VISCOSITY_FLUID/(PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID);
            GAMMA = 2 * ((HEAT_CAPACITY_FLUID - 1)/HEAT_CAPACITY_FLUID)*(1/H);
            
            D = [1 1 0 ; 1 1 0 ; 0 0 0];
            
            D_S = ((1 - POROSITY)^2/POROSITY)*D;
            D_SF = (1 - POROSITY)*D;
            D_F = POROSITY*D;
            
            stiffnessMatrix_S = sparse(8,8);
            stiffnessMatrix_SF = sparse(8,8);
            stiffnessMatrix_F = sparse(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, ~, Be, ~, J] = computeShapeFunction(obj, xi, eta);
                    stiffnessMatrix_S=stiffnessMatrix_S+(w(i)*w(j)*det(J)*transpose(Be)*(D_S*Be));
                    stiffnessMatrix_F=stiffnessMatrix_F+(w(i)*w(j)*det(J)*transpose(Be)*(D_F*Be));
                    stiffnessMatrix_SF=stiffnessMatrix_SF+(w(i)*w(j)*det(J)*transpose(Be)*(D_SF*Be));
                end
            end
            
            K2 = [stiffnessMatrix_S stiffnessMatrix_SF; stiffnessMatrix_SF stiffnessMatrix_F];
            
            ElementDampingMatrix=sparse(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(obj, xi, eta);
                    ElementDampingMatrix = ElementDampingMatrix + (w(i) * w(j) * transpose(N_mat) * N_mat * det(J));
                end
            end
            
            C = [ElementDampingMatrix, - ElementDampingMatrix; ...
                - ElementDampingMatrix, ElementDampingMatrix];
            
            dampingMatrix = beta * C + STANDARD_PRESSURE_FLUID * GAMMA * K2;
            
            dampingMatrix = obj.PorousStructureCoupling(dampingMatrix);
        end
        
         function dofs = getDofList(obj)
            dofs([1 3 5 7]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 4 6 8]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([9 11 13 15]) = obj.nodeArray.getDof('DISPLACEMENT_FLUID_X');
            dofs([10 12 14 16]) = obj.nodeArray.getDof('DISPLACEMENT_FLUID_Y');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,16); 
            vals([1 3 5 7]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6 8]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([9 11 13 15]) = obj.nodeArray.getDofValue('DISPLACEMENT_FLUID_X',step);
            vals([10 12 14 16]) = obj.nodeArray.getDofValue('DISPLACEMENT_FLUID_Y',step);
        end
        
        function vals = getFirstDerivativesVector(obj, step)
            vals = zeros(1,16);
            [~, vals([1 3 5 7]), ~] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, vals([9 11 13 15]), ~] = obj.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, vals([10 12 14 16]), ~] = obj.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(obj, step)
            vals = zeros(1,16);
            [~, ~, vals([1 3 5 7])] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, ~, vals([9 11 13 15])] = obj.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, ~, vals([10 12 14 16])] = obj.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
        end
        
        function update(obj)
        end
        
        function F = computeLocalForceVector(obj)
            F = sparse(1,16);
        end
        
        function matrix = PorousStructureCoupling(obj,matrix)
            for i=1:4
                check(i) = obj.nodeArray(i).getPropertyValue('IS_COUPLING_NODE');
            end
            
            if sum(check)==2
                couplingIds = find(check);
                couplingNodes = obj.nodeArray(couplingIds);
                matrixIds = (couplingIds -1)*2 +1;
                normalVector = obj.getNormalVector(couplingNodes);
                
                if abs(normalVector) == [1 0]
                    matrix(:,matrixIds) = matrix(:,matrixIds) + matrix(:,matrixIds+8);
                    matrix(matrixIds,:) = matrix(matrixIds,:) + matrix(matrixIds+8,:);
                    
                    matrix(:,matrixIds+8) = 0;
                    matrix(matrixIds+8,:) = 0;
                    
                elseif abs(normalVector) == [0 1]
                    matrix(:,matrixIds+1) = matrix(:,matrixIds+1) + matrix(:,matrixIds+9);
                    matrix(matrixIds+1,:) = matrix(matrixIds+1,:) + matrix(matrixIds+9,:);
                    
                    matrix(:,matrixIds+9) = 0;
                    matrix(matrixIds+9,:) = 0;
                    
                else
                    temp1 = matrix(:,matrixIds+8)*normalVector(1);
                    temp2 = matrix(:,matrixIds+9)*normalVector(2);
                    
                    matrix(:,matrixIds) = matrix(:,matrixIds) + temp1;
                    matrix(:,matrixIds+1) = matrix(:,matrixIds+1) + temp2;
                    
                    matrix(matrixIds,:) = matrix(matrixIds,:) + temp1.';
                    matrix(matrixIds+1,:) = matrix(matrixIds+1,:) + temp2.';
                    
                    matrix(:,matrixIds+8) = matrix(:,matrixIds+8) - temp1;
                    matrix(:,matrixIds+9) = matrix(:,matrixIds+9) - temp2;
                    
                    matrix(matrixIds+8,:) = matrix(matrixIds+8,:) - temp1.';
                    matrix(matrixIds+9,:) = matrix(matrixIds+9,:) - temp2.';
                    
                end
            elseif sum(check)==1 || sum(check)==3 || sum(check)==4
                e = MException('MATLAB:bm_mfem:cornerCoupling','Code can handle only 2-node coupling interfaces');
                throw(e);
            end
        end
        
    end
    
end
    

