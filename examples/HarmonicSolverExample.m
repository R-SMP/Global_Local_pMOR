%HARMONICSOLVEREXAMPLE An example how to perform a harmonic analysis in
% frequency range (frequency sweep analysis).
clear; close all

%create the model
[model, x0, xl, ~, ~] = createRectangularPlate(1, .2, 20, 8, 'elementType', 'ShellElement3d4n');
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', ...
    'ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
x0.fixAllDofs();

%add load condition
loadId = xl(end).getId();
model.addNewCondition('PointLoadCondition3d1n', 1, loadId);
model.getNode(loadId).setPropertyValue('POINT_LOAD',[0 0 -10]);

%define material properties
model.setPropertyValue('YOUNGS_MODULUS',2.1e11);
model.setPropertyValue('POISSON_RATIO',0.3);
model.setPropertyValue('DENSITY',7860);
model.setPropertyValue('THICKNESS',0.001);
model.setPropertyValue('NUMBER_GAUSS_POINT',2);

%create the solver
solver = HarmonicSolvingStrategy(model);

%define the frequency range in rad/s
freq = linspace(0.1,15,10);
for ii=freq
    solver.solve(ii);
end

%plot the results
result = model.getNode(loadId).getDofValue('DISPLACEMENT_Z','all');
semilogy(freq,abs(result(2:end)))