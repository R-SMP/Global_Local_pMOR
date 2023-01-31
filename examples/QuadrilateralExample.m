% A simple example with quadrilateral elements
clear
close all

% Dimension of the structure
Lx=10;
Ly=2;

% Number of elements in specific directions
nx=20;
ny=4;

% Model generation
[model, x0, xl, y0, yl] = createRectangularPlate(Lx, Ly, nx, ny, 'elementType', 'QuadrilateralElement2d4n');

% Assignment of material properties
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',30e6);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.0);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',3);
model.getAllElements.setPropertyValue('DENSITY',7860);

% Definition of Dofs and BCs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
x0.fixAllDofs();

% Definition of loading
addPointLoad(model.getNode((nx + 1)*(ny + 1)),300,[0 -1]);

% Solving
solver = SimpleSolvingStrategy(model);
x = solver.solve();

v = Visualization(model);
v.setScaling(1000);
v.plotUndeformed
v.plotDeformed
    
