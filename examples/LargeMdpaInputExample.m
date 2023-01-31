% LARGEMDPAINPUTEXAMPLE Read a large mdpa model and plot it
clear
close all

%specify input file
io = MdpaInput('large_mdpa_input.mdpa');
%read the model
model = io.readModel();

%set properties
model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
model.getAllElements.setPropertyValue('POISSON_RATIO',.3);
model.getAllElements.setPropertyValue('THICKNESS', 0.005);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',4);
model.getAllElements.setPropertyValue('DENSITY',7860);

%add dofs
model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
    "ROTATION_X", "ROTATION_Y",  "ROTATION_Z"]);

%set boundary conditions for all elements in 'left_support' and 'right_support':
model.getModelPart('GENERIC_left_support').getNodes().fixAllDofs();
model.getModelPart('GENERIC_right_support').getNodes().fixAllDofs();

%set load for all elements in 'inner_circle':
model.getModelPart('GENERIC_inner_circle').getNodes().setDofLoad('DISPLACEMENT_X',50000);
model.getModelPart('GENERIC_inner_circle').getNodes().setDofLoad('DISPLACEMENT_Y',50000);

%perform static solve
fprintf('Solving static analysis ... ');
tic
solver = SimpleSolvingStrategy(model);
solver.solve();
fprintf('done in %f seconds\n', toc);

%visualize
v=Visualization(model);
v.setScaling(1e2);
v.plotUndeformed()
v.plotDeformed
