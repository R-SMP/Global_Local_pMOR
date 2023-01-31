%SHELLEXAMPLE Example analyses for the shell element
%   A square plate meshed with 10 x 10 ShellElement3d4n is subjected to a
%   static loading on the middle node. The results can be visualized.
%
%   A eigenvalue analysis for the same plate is performed in the second
%   part. The resulting eigenmodes can be visualized using
%   v.plotDeformed(n), where n is the number of the mode to be visualized

%% Static analysis
close all;
clear; 

% % Initialization
nx = 10;
ny = 10;

fprintf('Creating a square plate with %u x %u elements\n',nx,ny);

[ model, x0, xl, y0, yl ] = createRectangularPlate( 1, 1, nx, ny,'elementType', 'ShellElement3d4n');
model.getAllNodes.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y','DISPLACEMENT_X', 'DISPLACEMENT_Y', 'ROTATION_Z'});

model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
model.getAllElements.setPropertyValue('POISSON_RATIO',.3);
model.getAllElements.setPropertyValue('THICKNESS', 0.005);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',4);
model.getAllElements.setPropertyValue('DENSITY',7860);

middle = fix((nx+1)*(ny+1)/2)+1; 

support = [x0 xl y0 yl];
support.fixAllDofs();

% % Static solve
model.getNode(middle).setDofLoad('DISPLACEMENT_Z',-1000);
fprintf('Solving static analysis ... ');
tic
solver = SimpleSolvingStrategy(model);
solver.solve();
fprintf('done in %f seconds\n', toc);
fprintf('Displacement at middle node: %f\n', ...
    model.getNode(middle).getDofValue('DISPLACEMENT_Z'))

% % Plot result
v=Visualization(model);
v.setScaling(1e2);
v.plotUndeformed()
v.plotDeformed()
view(-15,20)
zlim([-.4 .4])
title('Static displacement caused by load on the middle node')

str = input('Continue with eigenvalue analysis? Y/N [Y]: ','s');
if ~isempty(str) && ~strcmpi(str,'y')
    return
end

%% Eigenvalue analysis
close all;
clear; 

% % Initialization
nx = 10;
ny = 10;

fprintf('Creating a square plate with %u x %u elements\n',nx,ny);

[ model, x0, xl, y0, yl ] = createRectangularPlate( 1, 1, nx, ny,'elementType', 'ShellElement3d4n');
model.getAllNodes.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y','DISPLACEMENT_X', 'DISPLACEMENT_Y', 'ROTATION_Z'});

model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
model.getAllElements.setPropertyValue('POISSON_RATIO',.3);
model.getAllElements.setPropertyValue('THICKNESS', 0.005);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',4);
model.getAllElements.setPropertyValue('DENSITY',7860);

support = [x0 xl y0 yl];
support.fixAllDofs();

% % Solve for eigenvalues
fprintf('Solving eigenvalue analysis ... ');
tic
solver = EigensolverStrategy(model); 
solver.solve(10);
solver.assignModeShapes;
fprintf('done in %f seconds\n', toc);
ef = solver.getEigenfrequencies('Hz');
fprintf(['Eigenfrequencies of the plate:\n', regexprep(num2str(ef'),'\s+',' Hz\n'), ' Hz\n'])

% % Plot results
v=Visualization(model);
v.plotUndeformed()
v.plotDeformed(3) % change number for the other eigenmodes
view(-15,20)
zlim([-.4 .4])
title('Third eigenmode')
