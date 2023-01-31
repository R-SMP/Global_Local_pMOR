% ANSYSINPUTGIDOUTPUTEXAMPLE Read a large ANSYS model and export the results to GiD
clear
close all

%check, if the ansys path is specified correctly
%see utilities/getPaths_changeme.m for further information
if ~(exist('getPaths','file') == 2)
    error('Script getPaths not found. Please create it to use ANSYS import.')
end
if ~(exist(getPaths('ansys'),'file') == 2)
	error('ANSYS executable not found. Please check your getPaths script.')
end

%specify input file
io = AnsysInput('examples/AnsysInputGiDOutputExample.inp', getPaths('ansys'));

%read the model 
%all properties, dofs, boundary conditions, loads, etc. are read from the 
%ANSYS input file
fprintf('Starting ANSYS and reading model ... ');
tic
model = io.readModel();
fprintf('done in %f seconds\n', toc);

%perform modal analysis
fprintf('Solving modal analysis ... ');
tic
solver = EigensolverStrategy(model);
solver.solve(10);
solver.assignModeShapes();
fprintf('done in %f seconds\n', toc);

%export the results to GiD
ofile = 'examples/AnsysInputGiDOutputExample';
out = GiDOutput(ofile);
out.setPreference('writeModelParts',true);
out.setPreference('nodalResults',["DISPLACEMENT"]);
out.writeMesh(model);
out.writeResults(model);
fprintf(['Output written to ' ofile '.res\n']);
