%SOILSUPERELEMENTEXAMPLE Large example for the SoilSuperElement
%
% Soilsuperelement has the same number of nodes in x and y direction (Ny=Nx)
%
% Chosing of parameters Nx, dx and eventually damping:
% dx: depends on frequency f [Hz], there should be at least 6 or 8 nodes in each wave length (avoid of aliasing)
%     dx <= sqrt(G/density)/f/6        % at least 6 nodes
%                                      % G = abs( E*(1-2*damping*1i)./(2*(1+Poisson)) )
%
% Nx: absolute upper limit = 147 (size of stiffnessmatrix = ( 3*(Nx/2)^2 ) x ( 3*(Nx/2)^2 ) <= program capacity)
%     lower limt: depends on dx, recommended Nx >= 64 (avoid of leakage effect)
%     checking of leakage effect: plots of displacement along x (and y) axis, no leakage effect means that displacement decays to zero at both
%     ends of x axis
%     if there is leakage effect, dx and/or Nx and/or damping of soil material must be raised
%
% Some reference values: f-dx-Nx-damping: 25-0.5-100-0.025, 25-0.3-100-0.05, 30-0.5-100-0.025, 30-0.2-100-0.05, 
%                                         60-0.3-100-0.025, 60-0.2-100-0.05, 120-0.1-100-0.05

clear;

% Input parameters
Nx = 64;
dx = 0.3;

f = 25; % Input is frequency (f); not circular frequency (omega);

% soil material
h1 = 1;
E1 = 45e6;
Poisson1 = 0.332;
density1 = 1800;
damping1 = 0.05;

% Half space
E2 = 45e6;
Poisson2 = 0.332;
density2 = 1800;
damping2 = 0.05;

% Initialization
nodeArray = generateSoilSuperElement(dx,Nx);

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

ele01 = SoilSuperElement(1,nodeArray,Nx);

% elementArray = [ele01];

% properties of upper layer
ele01.setPropertyValue('THICKNESS',h1);                              
ele01.setPropertyValue('YOUNGS_MODULUS',E1);
ele01.setPropertyValue('POISSON_RATIO',Poisson1);
ele01.setPropertyValue('DENSITY',density1);
ele01.setPropertyValue('ELEMENTAL_DAMPING',damping1);

% properties of half space
ele01.setPropertyValue('YOUNGS_MODULUS_HS',E2);        
ele01.setPropertyValue('POISSON_RATIO_HS',Poisson2);
ele01.setPropertyValue('DENSITY_HS',density2);
ele01.setPropertyValue('ELEMENTAL_DAMPING_HS',damping2);

% discretization parameters
ele01.setPropertyValue('SAMPLE_NUMBER',Nx);          
ele01.setPropertyValue('SAMPLE_RATE',dx);
ele01.setPropertyValue('FREQUENCY',f);

% get NodeId for point load through coordinate (x,y)
j = getLocalNodeId(nodeArray,0,0);
addPointLoad(nodeArray(j),6,[1 0 0]);

model = FemModel(nodeArray, ele01);

assembling = SimpleAssembler(model);
% calculation of dynamic stiffness matrix
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);      

% calculation of dynamic displacement (no static solution !!!!)
solver = SimpleSolvingStrategy(model);
displacement = solver.solve();              

% export results to GiD
out = GiDOutput('examples/SoilsuperelementGiDOutput');
out.setPreference('writeModelParts',false);
out.writeMesh(model);
out.setPreference('nodalResults',["DISPLACEMENT"]);
out.writeResults(model);




