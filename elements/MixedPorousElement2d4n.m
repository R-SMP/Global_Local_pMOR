classdef MixedPorousElement2d4n < PorousElement2d4n
    %   MixedPorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here

    %   DOF: Solid-Displacement & Pore-Pressure (see: ATALLA 1998)
    
    properties (Access = private)
    end
    
    methods
        
        % constructor
        function obj = MixedPorousElement2d4n(id,nodeArray)
            requiredPropertyNames = cellstr(["DENSITY_SOLID", ...
                "LAMBDA_SOLID", "MUE_SOLID", "DAMPING_SOLID", ...
                "DENSITY_FLUID", "VISCOSITY_FLUID", ...
                "STANDARD_PRESSURE_FLUID", "HEAT_CAPACITY_FLUID", ...
                "PRANDTL_NUMBER_FLUID", ...
                "POROSITY", "TORTUOSITY", "FLOW_RESISTIVITY", ...
                "VISCOUS_LENGHT", "THERMAL_LENGTH", ...
                "FREQUENCY", "NUMBER_GAUSS_POINT"]);
            
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
                "ACOUSTIC_PRESSURE"]);
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

            TORTUOSITY = obj.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = obj.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = obj.getPropertyValue('VISCOUS_LENGHT');
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * omega * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2))^0.5;
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            
            % Determine Coefficients for EMat
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (HEAT_CAPACITY_FLUID - 1)*((1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5))^-1);
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            GAMMA = POROSITY * ((EQUIVALENT_COUPLING_DENSITY /EQUIVALENT_FLUID_DENSITY) - (Q/R)) ;
            
            % Matrix valid for isotropic material only
            D = [LAMBDA + 2 * MUE,LAMBDA,0; ...
                LAMBDA,LAMBDA + 2 * MUE,0; ...
                0,0,MUE];
            
            K_Matrix = sparse(8,8);
            C_Matrix = sparse(8,4);
            H_Matrix = sparse(4,4);
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                        [N_mat, ~, Be, B, J] = computeShapeFunction(obj, zeta, eta);
                        K_Matrix=K_Matrix+(w(i)*w(j)*det(J)*transpose(Be)*(D*Be));
                        C_Matrix=C_Matrix+(w(i)*w(j)*det(J)*GAMMA*transpose(N_mat) * B);
                        H_Matrix=H_Matrix+(w(i)*w(j)*det(J)*(POROSITY^2/EQUIVALENT_FLUID_DENSITY)*transpose(B) * B);
                end
            end
            
            stiffnessMatrix = [K_Matrix, - C_Matrix; zeros(4,8), H_Matrix];
        end
            
        function massMatrix = computeLocalMassMatrix(obj)
            % Read Material Properties
            DENSITY_FLUID = obj.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = obj.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = obj.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = obj.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = obj.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = obj.getPropertyValue('POROSITY');
            THERMAL_LENGTH = obj.getPropertyValue('THERMAL_LENGTH');
            
            omega = obj.getPropertyValue('FREQUENCY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            DENSITY_SOLID = obj.getPropertyValue('DENSITY_SOLID');

            TORTUOSITY = obj.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = obj.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = obj.getPropertyValue('VISCOUS_LENGHT');
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * omega *  TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2))^0.5;
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_DENSITY = EQUIVALENT_SOLID_DENSITY - EQUIVALENT_COUPLING_DENSITY^2/EQUIVALENT_FLUID_DENSITY; 
            
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (HEAT_CAPACITY_FLUID - 1)*((1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5))^-1);
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            GAMMA = POROSITY * ((EQUIVALENT_COUPLING_DENSITY /EQUIVALENT_FLUID_DENSITY) - (Q/R)) ;
            
            M_Matrix = sparse(8,8);
            C_Matrix = sparse(8,4);
            Q_Matrix = sparse(4,4);
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                        [N_mat, N, ~, B, J] = computeShapeFunction(obj, zeta, eta);
                        M_Matrix=M_Matrix+(w(i)*w(j)*det(J)*EQUIVALENT_DENSITY*transpose(N_mat) * N_mat);
                        C_Matrix=C_Matrix+(w(i)*w(j)*det(J)*GAMMA*transpose(N_mat) * B);
                        Q_Matrix=Q_Matrix+(w(i)*w(j)*det(J)*(POROSITY^2/R)*transpose(N)*(N));
                end
            end
            
            massMatrix = [M_Matrix, sparse(8,4); transpose(C_Matrix), Q_Matrix];   
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            dampingMatrix = sparse(12,12);
        end
        
        function dofs = getDofList(obj)
            dofs([1 3 5 7]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 4 6 8]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([9 10 11 12]) = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,16); 
            vals([1 3 5 7]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6 8]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([9 10 11 12]) = obj.nodeArray.getDofValue('ACOUSTIC_PRESSURE',step);
        end
        
        function vals = getFirstDerivativesVector(obj, step)
            vals = zeros(1,16);
            [~, vals([1 3 5 7]), ~] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, vals([9 10 11 12]), ~] = obj.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(obj, step)
            vals = zeros(1,16);
            [~, ~, vals([1 3 5 7])] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, ~, vals([9 10 11 12])] = obj.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function update(obj)
        end
        
        function F = computeLocalForceVector(obj)
            F = sparse(1,16);
        end
        
    end
    
end
    

