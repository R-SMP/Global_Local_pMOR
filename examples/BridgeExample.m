%BRIDGEEXAMPLE Static analysis of a bridge model meshed with bars
%   21 BarElement3d2n with different cross sections are used in the model.
%   The model is taken from Felippa's IFEM ch. 21

clear
close all

%build the model
model = FemModel();

%add nodes (id, x, y, z)
model.addNewNode(1,0,0,0);
model.addNewNode(2,10,5,0);
model.addNewNode(3,10,0,0);
model.addNewNode(4,20,8,0);
model.addNewNode(5,20,0,0);
model.addNewNode(6,30,9,0);
model.addNewNode(7,30,0,0);
model.addNewNode(8,40,8,0);
model.addNewNode(9,40,0,0);
model.addNewNode(10,50,5,0);
model.addNewNode(11,50,0,0);
model.addNewNode(12,60,0,0);

%add dofs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

%add elements (type, id, nodes)
model.addNewElement('BarElement3d2n',1,[1 3]);
model.addNewElement('BarElement3d2n',2,[3 5]);
model.addNewElement('BarElement3d2n',3,[5 7]);
model.addNewElement('BarElement3d2n',4,[7 9]);
model.addNewElement('BarElement3d2n',5,[9 11]);
model.addNewElement('BarElement3d2n',6,[11 12]);
model.addNewElement('BarElement3d2n',7,[1 2]);
model.addNewElement('BarElement3d2n',8,[2 4]);
model.addNewElement('BarElement3d2n',9,[4 6]);
model.addNewElement('BarElement3d2n',10,[6 8]);
model.addNewElement('BarElement3d2n',11,[8 10]);
model.addNewElement('BarElement3d2n',12,[10 12]);
model.addNewElement('BarElement3d2n',13,[2 3]);
model.addNewElement('BarElement3d2n',14,[4 5]);
model.addNewElement('BarElement3d2n',15,[6 7]);
model.addNewElement('BarElement3d2n',16,[8 9]);
model.addNewElement('BarElement3d2n',17,[10 11]);
model.addNewElement('BarElement3d2n',18,[2 5]);
model.addNewElement('BarElement3d2n',19,[4 7]);
model.addNewElement('BarElement3d2n',20,[7 8]);
model.addNewElement('BarElement3d2n',21,[9 10]);

%set properties
model.getElements(18:21).setPropertyValue('CROSS_SECTION',1);
model.getElements(1:6).setPropertyValue('CROSS_SECTION',2);
model.getElements(7:12).setPropertyValue('CROSS_SECTION',10);
model.getElements(13:17).setPropertyValue('CROSS_SECTION',3);
model.getAllElements().setPropertyValue('YOUNGS_MODULUS',1000);

%set boundary conditions
model.getNode(1).fixAllDofs();
model.getNode(12).fixDof('DISPLACEMENT_Y');
model.getAllNodes().fixDof('DISPLACEMENT_Z');

%set loads using function addPointLoad from utilities (nodes, magnitude,
%direction)
addPointLoad(model.getNodes([3 5 9 11]),10,[0 -1 0]);
addPointLoad(model.getNode(7),16,[0 -1 0]);

%perform the static analysis
solver = SimpleSolvingStrategy(model);
solver.solve();

%visualize the result
v=Visualization(model);
v.plotUndeformed()
v.plotDeformed()
