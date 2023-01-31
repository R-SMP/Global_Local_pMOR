%% COUPLING OF SOIL AND GABION

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERARTION OF SOILSUPERELEMENT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parameters
Nx = 64;
dx = 0.2;
f = 120; % Input is frequency (f); not circular frequency (omega);

% Initialization
nodeArray = generateSoilSuperElement(dx,Nx);
ele01 = SoilSuperElement(1,nodeArray,Nx);

% Properties of upper layer
ele01.setPropertyValue('THICKNESS',1);                              
ele01.setPropertyValue('YOUNGS_MODULUS',45e6);
ele01.setPropertyValue('POISSON_RATIO',0.332);
ele01.setPropertyValue('DENSITY',1800);
ele01.setPropertyValue('ELEMENTAL_DAMPING',0.05);

% Properties of half space
ele01.setPropertyValue('YOUNGS_MODULUS_HS',45e6);        
ele01.setPropertyValue('POISSON_RATIO_HS',0.332);
ele01.setPropertyValue('DENSITY_HS',1800);
ele01.setPropertyValue('ELEMENTAL_DAMPING_HS',0.05);

% Discretization parameters
ele01.setPropertyValue('SAMPLE_NUMBER',Nx);          
ele01.setPropertyValue('SAMPLE_RATE',dx);
ele01.setPropertyValue('FREQUENCY',f);

model = FemModel(nodeArray, ele01);

%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF GABION %
%%%%%%%%%%%%%%%%%%%%%%%%

% Input (dimensions and position of gabion; discretization in z-direction)
Lx = 2;
Ly = 7;
Lz = 0.4;
px = 2;
py = 1;
nz = 2;

% Secondary parameters (NO INPUT)
nx = Lx/dx;
ny = Ly/dx;
dy = dx;
dz = Lz/nz;

% Generation of additional nodes 
id = length(nodeArray);
for k=1:(ny + 1)
    for j=2:(nz + 1)
        for i=1:(nx + 1)
            id=id+1;
            model.addNewNode(id,(i-1)*dx+px,(k-1)*dy+py,(j-1)*dz);
        end
    end
end

% Determination of interface-nodes
interfaceNodesU = zeros((ny +1),(nx +1));
for i=1:(Ly/dx)+1
    for j=1:(Lx/dx)+1
        interfaceNodesU(i,j) = (py/dx)*Nx+(i-1)*Nx + px/dx + j;
    end
end
interfaceNodesO = zeros((ny +1),(nx +1));
for i=1:ny+1
    for j=1:nx+1
        interfaceNodesO(i,j) = length(nodeArray) + j + (i-1)*(nx+1)*(nz); 
    end
end

% Generation of elements (flying gabion and interface)
id = 1;
for k=1:ny
    for j=1:nz-1
        for i=1:nx
            id=id+1;
            a = i + (j-1)*(nx+1) + (k-1)*(nx+1)*(nz) + length(nodeArray);
            model.addNewElement('HexahedronElement3d8n',id,[a, a+1, a+1+(nx+1)*(nz), a+(nx+1)*(nz), a+(nx+1), a+1+(nx+1), a+1+(nx+1)*(nz)+(nx+1), a+(nx+1)*(nz)+(nx+1)]);
        end
    end
end
id=id+1;
for i=1:ny
    for j=1:nx
        model.addNewElement('HexahedronElement3d8n',id,[interfaceNodesU(i,j), interfaceNodesU(i,j+1), interfaceNodesU(i+1,j+1), interfaceNodesU(i+1,j), interfaceNodesO(i,j), interfaceNodesO(i,j+1), interfaceNodesO(i+1,j+1), interfaceNodesO(i+1,j)]);
        id= id + 1;
    end
end

% Properties of gabion
elementArray = model.getAllElements;
elementArray = elementArray(2:end);
elementArray.setPropertyValue('YOUNGS_MODULUS',32);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',3);
elementArray.setPropertyValue('DENSITY',7860);

% Assignment DOFs (for gabion and soil-surface)
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF LOADING %
%%%%%%%%%%%%%%%%%%%%%%%%%

% get NodeId for point load through coordinate (x,y)
j = getLocalNodeId(nodeArray,5,4.5);
addPointLoad(nodeArray(j),6,[1 0 1]);

%%%%%%%%%%%
% SOLVING %
%%%%%%%%%%%

solver = HarmonicSolvingStrategy(model);
solver.solve(f*2*pi);             

% export results to GiD
out = GiDOutput('soilsuperelementGiDOutput');
out.setPreference('writeModelParts',false);
out.writeMesh(model);
out.setPreference('nodalResults',["DISPLACEMENT"]);
out.writeResults(model);




