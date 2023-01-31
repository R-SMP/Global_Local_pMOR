clear

io=MdpaInput('BlochHomogeneousInclusion.mdpa'); %specify input file.
model = io.readModel(); %read the model
model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',70e9);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);

%% Up to here everything is done as it would be set up of Standard FEM

solver = BlochInverse1D(model);

numberOfPhases = 60;
numberOfBands = 8;

[phases,frequencies]=solver.solve(numberOfPhases,numberOfBands);

%% Visualize results
figure()

for i=1:numberOfBands
    plot(phases,frequencies(i,:),'r')
    hold on
end

xlabel('Phase [rad]')
ylabel('Frequency [Hz]')
ylim([0 1e3])
xlim([0 pi])