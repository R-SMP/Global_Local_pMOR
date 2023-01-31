classdef AcousticElement2d3n < TriangularElement
    %   AcousticElement2d4n A triangular acoustic element
    %   The physical properties of air are initialized by default. Setting
    %   the flag IS_PML to true converts the element to a PML. The distance
    %   from the physical domain and the normal direction can be computed
    %   in a preprocessing step by setting PML_PREPROCESS to true.
    
    %   For fluid implementation see: FRANCK 2008: Finite-Elemente-Methoden, 
    %       Loesungsalgorithmen und Werkzeuge fuer die akustische Simulationstechnik
    %   For locally defined PML see BERIOT 2016: Perfectly matched layers
    %       for time harmonic acoustics
    
    properties (Access = protected)
    end
    
    methods
        % constructor
        function obj = AcousticElement2d3n(id, nodeArray)
            requiredProperties = cellstr(["YOUNGS_MODULUS","DENSITY","NUMBER_GAUSS_POINT","IS_PML","PML_PREPROCESS"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 3 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredProperties};
            end
            
            % call the super class constructor
            obj@TriangularElement(super_args{:});
            obj.dofNames = cellstr("ACOUSTIC_PRESSURE");
        end
        
      
        
        %Initialization
        function initialize(obj)            
            % Default setting for parameters of fluid (AIR) and
            % NUMBER_GAUSS_POINT
            obj.eProperties.setValue('DENSITY',1.21);
            obj.eProperties.setValue('YOUNGS_MODULUS',142771);
            obj.eProperties.setValue('NUMBER_GAUSS_POINT',3);
        end
        
        % member functions
        function [N, B, J] = computeShapeFunction(obj,tCoord)
            % Shape Function and Derivatives
            N = tCoord; %[xi, eta, 1-xi-eta]
            N_Diff_Par = [1 0 -1; 0 1 -1];
            
            % Coordinates
            ele_coords = zeros(3,2);
            for i=1:3
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            if obj.eProperties.getValue('IS_PML') && (~obj.eProperties.getValue('PML_PREPROCESS'))
                for ii=1:3
                    pml_dist = obj.nodeArray(ii).getPropertyValue('PML_COORDINATE');
                    pml_dir = obj.nodeArray(ii).getPropertyValue('PML_DIRECTION');
                    ele_coords(ii,1) = ele_coords(ii,1) - (1/1i)*pml_dist(1)*pml_dir(1);
                    ele_coords(ii,2) = ele_coords(ii,2) - (1/1i)*pml_dist(2)*pml_dir(2);
                end
            end
            
            % Jacobian
            J = N_Diff_Par * ele_coords;
            
            % Calculation of B-Matrix
            B=J\N_Diff_Par;
            
            if obj.eProperties.getValue('PML_PREPROCESS')
                N_Diff = B;
                B = sparse(3,6);
                B(1,1:2:end) = N_Diff(1,:);
                B(3,2:2:end) = N_Diff(1,:);
                B(2,2:2:end) = N_Diff(2,:);
                B(3,1:2:end) = N_Diff(2,:);
            end
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            rho = obj.tryGetPropertyValue('DENSITY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            if obj.eProperties.getValue('PML_PREPROCESS')
                stiffnessMatrix = sparse(6,6);
                factor = 1.0;
            else
                stiffnessMatrix=sparse(3,3);
                factor = 1/rho;
            end
            
            for ii=1:p
                [w,g] = returnGaussPointTrig(p,ii);
                [~, B, J] = computeShapeFunction(obj,g);
                stiffnessMatrix = stiffnessMatrix + ...
                    factor * w * (B.' * B) * 0.5 * det(J);
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            if obj.eProperties.getValue('PML_PREPROCESS')
                massMatrix=sparse(6,6);
                return
            end
            
            massMatrix=sparse(3,3);
            E = obj.tryGetPropertyValue('YOUNGS_MODULUS');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            for ii=1:p
                [w,g] = returnGaussPointTrig(p,ii);
                [N, ~, J] = computeShapeFunction(obj,g);
                massMatrix = massMatrix + ...
                    (1/E) * w * (N.' * N) * 0.5 * det(J);
            end
        end

        function dofs = getDofList(obj)
            if obj.eProperties.getValue('PML_PREPROCESS')
                dofs([1 3 5]) = obj.nodeArray.getDof('PML_DISTANCE_X'); 
                dofs([2 4 6]) = obj.nodeArray.getDof('PML_DISTANCE_Y');
            else
                dofs([1 2 3]) = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
            end
        end
        
        function vals = getValuesVector(obj, step)
            if obj.eProperties.getValue('PML_PREPROCESS')
                vals = zeros(1,6);
                vals([1 3 5]) = obj.nodeArray.getDofValue('PML_DISTANCE_X',step);
                vals([2 4 6]) = obj.nodeArray.getDofValue('PML_DISTANCE_Y',step);
            else
                vals = zeros(1,3);
                vals([1 2 3]) = obj.nodeArray.getDofValue('ACOUSTIC_PRESSURE',step);
            end
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,3);
            [~, vals([1 2 3]), ~] = element.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,3);
            [~, ~, vals([1 2 3])] = element.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function update(obj)
        end
        
        function postprocess(obj)
            if ~obj.eProperties.getValue('PML_PREPROCESS')
                return
            end
            phi = obj.getValuesVector('end');            
            points = [1 0 0;0 1 0;0 0 1]; % evaluate at all Gauss points
            n = size(points,1);
            
            for j = 1:n
                [~, B, ~] = obj.computeShapeFunction(points(j,:));
                res = B * phi';
                res = reshape(res,1,[]);
                res(3) = 0; % this is a 2d element
                
                for ii=1:n
                    tmp = obj.nodeArray(ii).getPropertyValue('PML_DIRECTION');
                    tmp = tmp + res;
                    obj.nodeArray(ii).setPropertyValue('PML_DIRECTION',tmp);
                end
            end
            
%             for ii=1:n
%                 if ~any(obj.nodeArray(ii).getPropertyValue('PML_DIRECTION'))
%                     disp('yo')
%                 end
%             end
        end
        
    end
    
end


