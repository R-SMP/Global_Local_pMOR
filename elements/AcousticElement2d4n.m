classdef AcousticElement2d4n < Element
    %   AcousticElement2d4n A quadrilateral acoustic element
    %   The physical properties of air are initialized by default. Setting
    %   the flag IS_PML to true converts the element to a PML. The distance
    %   from the physical domain and the normal direction can be computed
    %   in a preprocessing step by setting PML_PREPROCESS to true.
    
    %   For fluid implementation see: FRANCK 2008: Finite-Elemente-Methoden, 
    %       Loesungsalgorithmen und Werkzeuge fuer die akustische Simulationstechnik
    %   For locally defined PML see BERIOT 2016: Perfectly matched layers
    %       for time harmonic acoustics
    
    properties (Access = protected)
        lengthX
        lengthY
        lenghtZ
    end
    
    methods
        % constructor
        function obj = AcousticElement2d4n(id, nodeArray)
            requiredProperties = cellstr(["YOUNGS_MODULUS","DENSITY","NUMBER_GAUSS_POINT","IS_PML","PML_PREPROCESS"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredProperties};
            end
            
            % call the super class constructor
            obj@Element(super_args{:});
            obj.dofNames = cellstr(["ACOUSTIC_PRESSURE"]); %,"PML_DISTANCE_X", "PML_DISTANCE_Y"
        end
        
        function c = barycenter(obj)
            %ACOUSTICELEMENT2D4N.BARYCENTER returns the barycenter of the
            %element in the x-y plane
            diag1X = [obj.nodeArray(1).getX() obj.nodeArray(3).getX()];
            diag1Y = [obj.nodeArray(1).getY() obj.nodeArray(3).getY()];
            diag2X = [obj.nodeArray(2).getX() obj.nodeArray(4).getX()];
            diag2Y = [obj.nodeArray(2).getY() obj.nodeArray(4).getY()];
            
            try
                [c(1),c(2)] = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            catch e
                if (strcmp(e.identifier,'MATLAB:subsassigndimmismatch'))
                    msg = ['Barycenter of element ', ...
                        num2str(obj.getId()), ...
                        ' could not be computed. ', ...
                        'Check the element for convexity.'];
                    causeException = MException('MATLAB:bm_mfem:barycenter',msg);
                    e = addCause(e,causeException);
                elseif (strcmp(e.identifier,'MATLAB:UndefinedFunction'))
                    msg = ['Function ''polyxpoly'' not found. ', ...
                        'Please install the Mapping Toolbox.'];
                    causeException = MException('MATLAB:bm_mfem:toolboxnotfound',msg);
                    e = addCause(e,causeException);
                end
                rethrow(e);
            end
        end
        
        function c = checkConvexity(obj)
            %ACOUSTICELEMENT2D4N.CHECKCONVEXITY returns true, if the
            %element is convex
            diag1X = [obj.nodeArray(1).getX() obj.nodeArray(3).getX()];
            diag1Y = [obj.nodeArray(1).getY() obj.nodeArray(3).getY()];
            diag2X = [obj.nodeArray(2).getX() obj.nodeArray(4).getX()];
            diag2Y = [obj.nodeArray(2).getY() obj.nodeArray(4).getY()];
            
            c = 0;
            try
                res = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            catch e
                if (strcmp(e.identifier,'MATLAB:subsassigndimmismatch'))
                    c = 0;
                elseif (strcmp(e.identifier,'MATLAB:UndefinedFunction'))
                    msg = ['Function ''polyxpoly'' not found. ', ...
                        'Please install the Mapping Toolbox.'];
                    causeException = MException('MATLAB:bm_mfem:toolboxnotfound',msg);
                    e = addCause(e,causeException);
                    rethrow(e);
                end
                
            end
            
            if ~isempty(res); c = 1; end
        end
        
        %Initialization
        function initialize(obj)
            obj.lengthX = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(2).getCoords);
            
            obj.lengthY = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(4).getCoords);
            
            checkConvexity(obj);
            
            % Default setting for parameters of fluid (AIR) and
            % NUMBER_GAUSS_POINT
            obj.eProperties.setValue('DENSITY',1.21);
            obj.eProperties.setValue('YOUNGS_MODULUS',142771);
            obj.eProperties.setValue('NUMBER_GAUSS_POINT',2);
        end
        
        % member functions
        function [N, B, J] = computeShapeFunction(obj,xi,eta)
            % Shape Function and Derivatives
            N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4;
                -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];

            % Coordinates of the nodes forming one element
            ele_coords = zeros(4,2);
            for ii=1:4
                ele_coords(ii,1) = obj.nodeArray(ii).getX;
                ele_coords(ii,2) = obj.nodeArray(ii).getY;
            end
            
            if obj.eProperties.getValue('IS_PML') && (~obj.eProperties.getValue('PML_PREPROCESS'))
                for ii=1:4
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
                B = sparse(3,8);
                B(1,1:2:end) = N_Diff(1,:);
                B(3,2:2:end) = N_Diff(1,:);
                B(2,2:2:end) = N_Diff(2,:);
                B(3,1:2:end) = N_Diff(2,:);
            end
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            rho = obj.tryGetPropertyValue('DENSITY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            if obj.eProperties.getValue('PML_PREPROCESS')
                stiffnessMatrix = sparse(8,8);
                factor = 1.0;
            else
                stiffnessMatrix=sparse(4,4);
                factor = 1/rho;
            end
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, B, J] = computeShapeFunction(obj, xi, eta);
                    stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*factor*det(J)*transpose(B)*B);
                end
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            if obj.eProperties.getValue('PML_PREPROCESS')
                massMatrix=sparse(8,8);
                return
            end
            massMatrix=sparse(4,4);
            
            E = obj.tryGetPropertyValue('YOUNGS_MODULUS');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N, ~, J] = computeShapeFunction(obj, xi, eta);
                    massMatrix=massMatrix+(w(i)*w(j)*(1/E)*transpose(N)*N*det(J));
                end
            end
        end

        function dofs = getDofList(obj)
            if obj.eProperties.getValue('PML_PREPROCESS')
                dofs([1 3 5 7]) = obj.nodeArray.getDof('PML_DISTANCE_X'); 
                dofs([2 4 6 8]) = obj.nodeArray.getDof('PML_DISTANCE_Y');
            else
                dofs([1 2 3 4]) = obj.nodeArray.getDof('ACOUSTIC_PRESSURE');
            end
        end
        
        function vals = getValuesVector(obj, step)
            if obj.eProperties.getValue('PML_PREPROCESS')
                vals = zeros(1,8);
                vals([1 3 5 7]) = obj.nodeArray.getDofValue('PML_DISTANCE_X',step);
                vals([2 4 6 8]) = obj.nodeArray.getDofValue('PML_DISTANCE_Y',step);
            else
                vals = zeros(1,4);
                vals([1 2 3 4]) = obj.nodeArray.getDofValue('ACOUSTIC_PRESSURE',step);
            end
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,4);
            [~, vals([1 2 3 4]), ~] = element.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,4);
            [~, ~, vals([1 2 3 4])] = element.nodeArray.getDof('ACOUSTIC_PRESSURE').getAllValues(step);
        end
        
        function pl = draw(obj)
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(4).getX,...
                obj.nodeArray(1).getX];
            
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(4).getY, ...
                obj.nodeArray(1).getY];
            
            if(all(obj.getNodes().getDimension == 3))
                z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                    obj.nodeArray(3).getZ, obj.nodeArray(4).getZ, ...
                    obj.nodeArray(1).getZ];
                
                pl = line(x,y,z);
            else
                pl = line(x,y);
            end
        end
        
        function update(obj)
        end
        
        function postprocess(obj)
            if ~obj.eProperties.getValue('PML_PREPROCESS')
                return
            end
            phi = obj.getValuesVector('end');            
            points = [-1 -1;1 -1;1 1;-1 1]; % evaluate at all Gauss points
            n = size(points,1);
            
            for j = 1:n
                [~, B, ~] = obj.computeShapeFunction(points(j,1),points(j,2));
                res = B * phi';
                res = reshape(res,1,[]);
                res(3) = 0; % this is a 2d element
                
                for ii=1:n
                    tmp = obj.nodeArray(ii).getPropertyValue('PML_DIRECTION');
                    tmp = tmp + res;
                    obj.nodeArray(ii).setPropertyValue('PML_DIRECTION',tmp);
                end
            end
        end
        
    end
    
    % Definig element type for GIDOutput
    methods (Static)
        function o = getElementType()
            o = 'Quadrilateral';
        end
    end
    
end


