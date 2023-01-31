%BlochInverse1D    
%Unit-cell created and meshed with GiD
%%%% unit cell with l = 0.2m, h = 0.15m; a round inclusion in the center
%%% with r = 0.025m; a horizontal spring mass system will be added in the
%%% inclusion
clear;

io=MdpaInput('tests/input_data/Federtests_Kreis.mdpa'); %specify input file   
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

%Element properties for homogenous unic-cell: here aluminium
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);


%% Adding a mass-spring-system 
% with two springs and a mass


allNodes = model.getAllNodes();
massNodeID = length(allNodes)+1;

allElements = model.getAllElements();
massElementID = length(allElements)+1;
spring1ElementID = massElementID+1;
spring2ElementID = massElementID+2;


springNodes = model.getModelPart('GENERIC_0Grad').getNodes(); %Pre-defined nodes for the springs at 
%'GENERIC22,5Grad', 'GENERIC_45Grad', 'GENERIC_67.5Grad', 'GENERIC_90Grad'
%would also work for this example

leftSpringNode = springNodes(1,1);
leftSNId = getId(leftSpringNode);

rightSpringNode = springNodes(1,2);
rightSNId = getId(rightSpringNode);


% %% 2 springs, 1 mass
model.addNewNode(massNodeID,0.1,0.075,0); %Here: Center of the Unit-Cell
model.addNewElement('SpringDamperElement3d2n',spring1ElementID,[leftSNId massNodeID]);
model.addNewElement('SpringDamperElement3d2n',spring2ElementID,[massNodeID rightSNId]);
model.addNewElement('ConcentratedMassElement3d1n',massElementID, massNodeID);


model.getNode(leftSNId).addDof("DISPLACEMENT_Z"); %'SpringDamperElement3d2n' requieres DOF in z-direction (as well as x & y)
model.getNode(rightSNId).addDof("DISPLACEMENT_Z"); %'SpringDamperElement3d2n' requieres DOF in z-direction (as well as x & y)
model.getNode(massNodeID).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]); %'ConcentratedMassElement3d1n' requieres DOF in x,y,z-direction


model.getNode(leftSNId).fixDof('DISPLACEMENT_Z'); 
model.getNode(rightSNId).fixDof('DISPLACEMENT_Z');
model.getNode(massNodeID).fixDof('DISPLACEMENT_Z'); %z-displacements need to be fixed
model.getNode(massNodeID).fixDof('DISPLACEMENT_Y'); %has to be removed for 'GENERIC22,5Grad', 'GENERIC_45Grad', 'GENERIC_67.5Grad' and 'GENERIC_90Grad' 
%if the upper line (fixDOF_Y) is removed, the rigid body mode can be seen

model.getElement(spring1ElementID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8); %stiffness of the spring 1 
model.getElement(spring1ElementID).setPropertyValue('ELEMENTAL_DAMPING',0); %damping is not considered
model.getElement(spring2ElementID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);  %stiffness of the spring 2
model.getElement(spring2ElementID).setPropertyValue('ELEMENTAL_DAMPING',0);

model.getElement(massElementID).setPropertyValue('ELEMENTAL_MASS',7e0);
%% Fixing displacement perpendicular to the springs

%%% if this is not done, the rigid body mode can be seen

%%% this can also be done for horizontal springs with  
%%% model.getNode(massNodeID).fixDof('DISPLACEMENT_Y'); see above
%%% and for vertical springs
%%% model.getNode(massNodeID).fixDof('DISPLACEMENT_X');

% fixingNode = model.getModelPart('GENERIC_fixedNodes').getNodes();
% %'GENERIC_fixedNodes' are nodes, that are used as supports for one side
% %of the springs
% fixingNodeID = getId(fixingNode);
% fixingSpringID = massElementID + 3;
% 
% model.getNode(fixingNodeID).addDof("DISPLACEMENT_Z");
% 
% model.getNode(fixingNodeID).fixDof('DISPLACEMENT_X');
% model.getNode(fixingNodeID).fixDof('DISPLACEMENT_Y');
% model.getNode(fixingNodeID).fixDof('DISPLACEMENT_Z');
% 
% model.addNewElement('SpringDamperElement3d2n',fixingSpringID,[fixingNodeID massNodeID]);
% 
% model.getElement(fixingSpringID).setPropertyValue('ELEMENTAL_STIFFNESS',10e20);
% %the stiffness of the fixing springs has to be much higher than the
% stiffness of the springs belonging to the spring mass system
% model.getElement(fixingSpringID).setPropertyValue('ELEMENTAL_DAMPING',0);

%% visualization and assembling

v=Visualization(model); %set up visualization
v.plotUndeformed()  %visualize

%% Up to here everything is done as it would be set up of Standard FEM

% Create Bloch Solver
solver = BlochDirect1DSolvingStrategy(model);

% call the solve function of the solver
% maxFreq: maximal frequency 
% NumberOfFrequencies: how many frequency steps from 0 to max frequency
frequencies=linspace(0,5000,200);
[eigenValues]=solver.solve(frequencies);

%% Visualize results
figure()
for i=1:size(eigenValues,2)
scatter(real(eigenValues(:,i)),imag(eigenValues(:,i)))
end
title('Wave solutions')
xlabel('Real part')
ylabel('Imaginary part')
xlim([-2 2])
ylim([-2 2])


figure()
for i=1:size(eigenValues,2)
    freqs=frequencies(i)*ones(size(eigenValues(:,i)));
    phase=atan(imag(eigenValues(:,i))./real(eigenValues(:,i)));
    decay=log(abs(eigenValues(:,i)));
    
    in1=abs(eigenValues(:,i))<10;
    in2=abs(eigenValues(:,i))>0.1;
    in=find(and(in1,in2));
    decay=decay(in);
    phase=phase(in);
    freqs=freqs(in);
%     decay=decay(-10<(abs(eigenValues(:,i))-1));
%     phase=phase(-10<(abs(eigenValues(:,i))-1)> -10);
%     freqs=freqs(-10<(abs(eigenValues(:,i))-1)> -10);
    subplot(2,1,1)
    scatter(phase,freqs)
    xlabel('rad/m')
    hold on
    subplot(2,1,2)
    scatter(decay,freqs)
    xlabel('rad/m')
    hold on
end






