classdef ClassicalPorousElement2d4n < PorousElement2d4n
    %   ClassicalPorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Fluid-Displacement (see: KANG 1995 / RUMPLER 2012)
    
    %   Coupling between a ClassicalPorousElement and a structure is
    %   considered here, but yet works only for 2-node interfaces between elements 
    %   belonging to different media (not 3 or 1 node interfaces)

    
    properties (Access = private)
    end
    
    methods
        
        % constructor
        function obj = ClassicalPorousElement2d4n(id,nodeArray)
            requiredPropertyNames = cellstr(["DENSITY_SOLID", "LAMBDA_SOLID", ...
                "MUE_SOLID", "DAMPING_SOLID", "DENSITY_FLUID", ...
                "VISCOSITY_FLUID", "STANDARD_PRESSURE_FLUID", ...
                "HEAT_CAPACITY_FLUID", "PRANDTL_NUMBER_FLUID", "POROSITY", ...
                "TORTUOSITY", "FLOW_RESISTIVITY", "VISCOUS_LENGHT", ...
                "THERMAL_LENGTH", "FREQUENCY", "NUMBER_GAUSS_POINT"]);
            
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
            
            DENSITY_FLUID = obj.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = obj.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = obj.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = obj.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = obj.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = obj.getPropertyValue('POROSITY');
            THERMAL_LENGTH = obj.getPropertyValue('THERMAL_LENGTH');
            
            omega = obj.getPropertyValue('FREQUENCY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine Coefficients for EMat
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (HEAT_CAPACITY_FLUID - 1)*(1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5)^-1);
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            A = LAMBDA + ((1 - POROSITY)^2/ POROSITY) * K;
            
            % Matrix valid for isotropic material only
            D_S = [A + 2 * MUE,A,0; ...
                A,A + 2 * MUE,0; ...
                0,0,MUE];
            D_SF = [Q,Q,0; ...
                Q,Q,0; ...
                0,0,0];
            D_F = [R,R,0; ...
                R,R,0; ...
                0,0,0];
            
            stiffnessMatrix_S = sparse(8,8);
            stiffnessMatrix_F = sparse(8,8);
            stiffnessMatrix_SF = sparse(8,8);
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
            
            omega = obj.getPropertyValue('FREQUENCY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * sqrt(1 + ...
                4 * 1i *omega* TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2));
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            
            ElementmassMatrix = sparse(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(obj, xi, eta);
                    ElementmassMatrix = ElementmassMatrix + (w(i) * w(j) * transpose(N_mat) * N_mat * det(J));
                end
            end
            
            massMatrix = [EQUIVALENT_SOLID_DENSITY * ElementmassMatrix, EQUIVALENT_COUPLING_DENSITY * ElementmassMatrix; ...
                EQUIVALENT_COUPLING_DENSITY * ElementmassMatrix, EQUIVALENT_FLUID_DENSITY * ElementmassMatrix];
            
            for i=1:4
                check(i) = obj.nodeArray(i).getPropertyValue('IS_COUPLING_NODE');
            end
            
            if sum(check)==2
                couplingIds = find(check);
                couplingNodes = obj.nodeArray(couplingIds);
                matrixIds = (couplingIds -1)*2 +1;
                normalVector = obj.getNormalVector(couplingNodes);
                
                if abs(normalVector) == [1 0]
                    massMatrix(:,matrixIds) = massMatrix(:,matrixIds) + massMatrix(:,matrixIds+8);
                    massMatrix(matrixIds,:) = massMatrix(matrixIds,:) + massMatrix(matrixIds+8,:);
                    
                    massMatrix(:,matrixIds+8) = 0;
                    massMatrix(matrixIds+8,:) = 0;
                    
                elseif abs(normalVector) == [0 1]
                    massMatrix(:,matrixIds+1) = massMatrix(:,matrixIds+1) + massMatrix(:,matrixIds+9);
                    massMatrix(matrixIds+1,:) = massMatrix(matrixIds+1,:) + massMatrix(matrixIds+9,:);
                    
                    massMatrix(:,matrixIds+9) = 0;
                    massMatrix(matrixIds+9,:) = 0;
                    
                else
                    temp1 = massMatrix(:,matrixIds+8)*normalVector(1);
                    temp2 = massMatrix(:,matrixIds+9)*normalVector(2);
                    
                    massMatrix(:,matrixIds) = massMatrix(:,matrixIds) + temp1;
                    massMatrix(:,matrixIds+1) = massMatrix(:,matrixIds+1) + temp2;
                    
                    massMatrix(matrixIds,:) = massMatrix(matrixIds,:) + temp1.';
                    massMatrix(matrixIds+1,:) = massMatrix(matrixIds+1,:) + temp2.';
                    
                    massMatrix(:,matrixIds+8) = massMatrix(:,matrixIds+8) - temp1;
                    massMatrix(:,matrixIds+9) = massMatrix(:,matrixIds+9) - temp2;
                    
                    massMatrix(matrixIds+8,:) = massMatrix(matrixIds+8,:) - temp1.';
                    massMatrix(matrixIds+9,:) = massMatrix(matrixIds+9,:) - temp2.';
                    
                end
            elseif sum(check)==1 || sum(check)==3 || sum(check)==4
                e = MException('MATLAB:bm_mfem:cornerCoupling','Code can handle only 2-node coupling interfaces');
                throw(e);
            end
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            dampingMatrix = sparse(16,16);
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
        
    end
    
end
    

