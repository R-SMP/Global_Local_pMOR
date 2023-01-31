% A simple example which demonstrates how the coupling between a
% ClassicalPorousElement2d4n and a structure works

% Coupling yet works only for 2-node interfaces between elements belonging
% to different media (not 3 or 1 node interfaces)

clear
close all

% Dimension of the structure
node01 = Node(1,0,0,0);
node02 = Node(2,1,0,0);
node03 = Node(3,1,1,0);
node04 = Node(4,0,1,0);
node05 = Node(5,2,0,0);
node06 = Node(6,2,1,0);

nodeArray = [node01 node02 node03 node04 node05 node06];

nodeArray(1:6).addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
nodeArray(5:6).addDof({'DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y'});
nodeArray(2:3).addDof({'DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y'});

ele01 = QuadrilateralElement2d4n(1,[node01 node02 node03 node04]);
ele02 = ClassicalPorousElement2d4n(2,[node02 node05 node06 node03]);

elementArray = [ele01 ele02];

node02.setPropertyValue('IS_COUPLING_NODE',1);
node03.setPropertyValue('IS_COUPLING_NODE',1);

elementArray.setPropertyValue('YOUNGS_MODULUS',32);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',3);
elementArray.setPropertyValue('DENSITY',7860);

elementArray.setPropertyValue('DENSITY_SOLID',750);
elementArray.setPropertyValue('LAMBDA_SOLID',487500);
elementArray.setPropertyValue('MUE_SOLID',325000);
elementArray.setPropertyValue('DAMPING_SOLID',0.0);

elementArray.setPropertyValue('POROSITY',0.96);
elementArray.setPropertyValue('TORTUOSITY',1.7);
elementArray.setPropertyValue('FLOW_RESISTIVITY',32e3);
elementArray.setPropertyValue('VISCOUS_LENGHT',90e-6);
elementArray.setPropertyValue('THERMAL_LENGTH',165e-6);

elementArray.setPropertyValue('FREQUENCY',200*2*pi);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
massMatrix = assembling.assembleGlobalMassMatrix(model);
