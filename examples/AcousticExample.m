% Acoustic Example
% Rectangular domain; rigid boundary on the left side; free radiation on
% all other boundaries;

clear

io = MdpaInput('AcousticExample.mdpa');
model = io.readModel;

omega = 2*pi*40;
% Definition of Properties of acoustic elements (undefined properties correspond to air by default)
model.getAllNodes.addDof({'ACOUSTIC_PRESSURE'});

robin_nodes = model.getModelPart('LinePressure2D_Robin').getNodes.getId;
robin_condition = model.getModelPart('LinePressure2D_Robin').getConditions.getId;
neumann_nodes = model.getModelPart('LinePressure2D_Neumann').getNodes.getId;
neumann_condition = model.getModelPart('LinePressure2D_Neumann').getConditions.getId;

%-DIRICHLET-
model.getNode(243).prescribeDof('ACOUSTIC_PRESSURE',1);

%-NEUMANN-
for i=1:length(neumann_condition)
    model.getCondition(neumann_condition(i)).setPropertyValue('FREQUENCY',omega);
end
for i=1:length(neumann_nodes)
    model.getNode(neumann_nodes(i)).setPropertyValue('NORMAL_DISPLACEMENT',0);
end

%-ROBIN-
for i=1:length(robin_nodes)
    model.getNode(robin_nodes(i)).setPropertyValue('ADMITTANCE',2.41e-3);
end

% Solving for unprescribed DOFS
solver = HarmonicSolvingStrategy(model);
solver.solve(omega);

% Generation of GID-output 
%export the results to GiD
out = GiDOutput('examples/AcousticExample');
out.setPreference('writeModelParts',false);
out.writeMesh(model);
out.setPreference('nodalResults',["ACOUSTIC_PRESSURE"]);
out.writeResults(model);
    
