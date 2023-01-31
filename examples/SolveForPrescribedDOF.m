% Solving for remaining DOF-values in case of given DOF-values
% One can only put loads on unprescribed DOFS!!!!

clear
close all

% Dimension of the structure
Lx=1;
Ly=1;

% Number of elements in specific directions
nx=5;
ny=5;

[model, x0, ~, ~, ~] = createRectangularPlate(Lx, Ly, nx, ny, 'elementType', 'QuadrilateralElement2d4n');

% Assignment of DOFs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

% assignment of material properties
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',30e6);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.0);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',3);
model.getAllElements.setPropertyValue('DENSITY',7860);

% Definition of BCs
x0.fixAllDofs();

% Prescription of DOFS
model.getNode((nx+1)*(ny+1)).prescribeDof('DISPLACEMENT_Y',-1.55555555555556e-07);

% % Checking result of prescription
% getValue(getDofArray(model.getNode(nx+1)))
  
% % % Definition of loading (ONLY Loads on unprescribed DOFS are possible!!!)
% addPointLoad(model.getNode(nx+1),1,[0 -1]);

% % Determination of global matrices
% assembling = SimpleAssembler(model);
% [stiffnessMatrix, Kred] = assembling.assembleGlobalStiffnessMatrix(model);           
% [massMatrix, Mred] = assembling.assembleGlobalMassMatrix(model);

% Solving for unprescribed DOFS
solver = PrescribedDofSolvingStrategy(model,0);
solver.solve();
displacement = model.getDofArray.getValue();
nodalForces = solver.getNodalForces();

% % Use this solver instead to check on results
% solver = SimpleSolvingStrategy(model);
% x = solver.solve();
% step = 1;
% VerschiebungDofs = model.getDofArray.getValue(step);
% nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(100000);
v.plotUndeformed
v.plotDeformed
    