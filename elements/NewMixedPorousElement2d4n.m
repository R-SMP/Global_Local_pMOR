classdef NewMixedPorousElement2d4n < PorousElement2d4n
    %   NewMixedPorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Pore-Pressure (see: ATALLA 2001)
    
    properties (Access = private)
    end
    
    methods
        
        % constructor
        function obj = NewMixedPorousElement2d4n(id,nodeArray)
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
            DENSITY_SOLID = obj.getPropertyValue('DENSITY_SOLID');
            
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
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_DENSITY = EQUIVALENT_SOLID_DENSITY - EQUIVALENT_COUPLING_DENSITY^2/EQUIVALENT_FLUID_DENSITY; 
            
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
            M_Matrix = sparse(8,8);
            H_Matrix = sparse(4,4);
            Q_Matrix = sparse(4,4);
            C1_Matrix = sparse(8,4);
            C2_Matrix = sparse(8,4);
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                        [N_mat, N, Be, B, J] = computeShapeFunction(obj, zeta, eta);
                        Bu = [B(1,1),B(2,1),B(1,2),B(2,2),B(1,3),B(2,3),B(1,4),B(2,4)];
                        K_Matrix=K_Matrix+(w(i)*w(j)*det(J)*transpose(Be)*(D*Be));
                        M_Matrix=M_Matrix+(w(i)*w(j)*det(J)*EQUIVALENT_DENSITY*transpose(N_mat) * N_mat);
                        H_Matrix=H_Matrix+(w(i)*w(j)*det(J)*(POROSITY^2/EQUIVALENT_FLUID_DENSITY)*transpose(B) * B);
                        Q_Matrix=Q_Matrix+(w(i)*w(j)*det(J)*(POROSITY^2/R)*transpose(N)*(N));
                        C1_Matrix=C1_Matrix+(w(i)*w(j)*det(J)*(GAMMA+1)*transpose(N_mat) * B);
                        C2_Matrix=C2_Matrix+(w(i)*w(j)*det(J)*transpose(Bu) * N);
%                         Using Approximation for Q and R (RUMPLER 2012)
%                         C1_Matrix=C1_Matrix+(w(i)*w(j)*det(J)*(GAMMA+POROSITY*(1+Q/R))*transpose(N_mat) * B);
%                         C2_Matrix=C2_Matrix+(w(i)*w(j)*det(J)*(POROSITY*(1+Q/R))*transpose(Bu) * N);
                end
            end
            stiffnessMatrix = [K_Matrix - omega^2 * M_Matrix ,- (C1_Matrix + C2_Matrix); - (omega^2)*transpose(C1_Matrix + C2_Matrix), H_Matrix - (omega^2)*Q_Matrix];
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            % massMatrix is integrated into dynamic stiffnessMatrix
            massMatrix = sparse(12,12);
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
    

