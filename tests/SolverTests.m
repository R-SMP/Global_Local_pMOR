classdef SolverTests <  matlab.unittest.TestCase
    %SOLVERTESTS Class for solver tests
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        function eigensolverTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %parameters
            mass = 1.0;
            stiffness = 10.0;
            
            %model
            model = FemModel();
            
            n01 = model.addNewNode(1,1,0,0);
            n02 = model.addNewNode(2,2,0,0);
            base = model.addNewNode(3,0,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            s01 = model.addNewElement('SpringDamperElement3d2n', 1, [3 1]);
            s01.setPropertyValue('ELEMENTAL_STIFFNESS', 2*stiffness);
            m01 = model.addNewElement('ConcentratedMassElement3d1n', 2, 1);
            m01.setPropertyValue('ELEMENTAL_MASS', 2*mass);
            s02 = model.addNewElement('SpringDamperElement3d2n', 3, [1 2]);
            s02.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            m02 = model.addNewElement('ConcentratedMassElement3d1n', 4, 2);
            m02.setPropertyValue('ELEMENTAL_MASS', mass);
            
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            base.fixDof('DISPLACEMENT_X');
            
            solver = EigensolverStrategy(model);
            solver.solve(2);
            solver.assignModeShapes();
            
            %eigenfrequencies in rad/s
            expectedEigenfrequencies = [sqrt(stiffness/(2*mass)) sqrt(2*stiffness/mass)]';
            actualEigenfrequencies = solver.getEigenfrequencies;
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                    'Within', AbsoluteTolerance(1e-7)))
        end
        
        function dampedHarmonicTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            a = .01;
            m = createBeam(.4,10,'IY',a^4/12,'IZ',a^4/12,'IT',a^4/12,'CROSS_SECTION',a^2);
            m.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z','ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
            
            m.getNode(1).fixAllDofs;
            m.getNode(11).setDofLoad('DISPLACEMENT_Z',100);
            m.getAllElements.addProperty('RAYLEIGH_ALPHA',4.64);
            m.getAllElements.addProperty('RAYLEIGH_BETA',9.1e-5);
            
            s = EigensolverStrategy(m);
            s.solveModalSuperposition(2*pi*logspace(1,4,100),10);
            
            disp = abs(m.getNode(11).getDofValue('DISPLACEMENT_Z','all'));
            load('tests/input_data/test_data.mat','damped_harmonic_disp');
            
            testCase.assertThat(disp, IsEqualTo(damped_harmonic_disp, ...
                'Within', AbsoluteTolerance(1e-5)))
            
        end
        
        function newmarkTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %parameters
            mass = 1.0;
            stiffness = 10.0;
            damping = 2.0;
            uInit = 0.1;
            vInit = 0.0;
            
            %analytical solution
            omega = sqrt(stiffness/mass);
            D = damping / (2 * mass * omega);
            omega_D = omega * sqrt(1 - power(D, 2));
            delta = damping / (2 * mass);
            theta = atan(- (vInit + uInit * delta) / (omega_D * uInit));
            A = sqrt(power(uInit, 2) + power((vInit + uInit * delta) / omega_D, 2));
            
            %model
            model = FemModel();
            
            n01 = model.addNewNode(1,0,1,0);
            n02 = model.addNewNode(2,0,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            springEle = model.addNewElement('SpringDamperElement3d2n', 1, [1 2]);
            springEle.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            springEle.setPropertyValue('ELEMENTAL_DAMPING', damping);
            massEle = model.addNewElement('ConcentratedMassElement3d1n', 2, 1);
            massEle.setPropertyValue('ELEMENTAL_MASS', mass);
            
            model.getAllNodes.fixDof('DISPLACEMENT_X');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            n02.fixDof('DISPLACEMENT_Y');
            
            n01.setDofValue('DISPLACEMENT_Y', uInit);
            
            dt = 0.1;
            time = 0;
            endTime = 4;
            solver = NewmarkSolvingStrategy(model, dt);
            while time < endTime
                solver.solve();
                time = time + dt;
                
                actualDisplacementY = n01.getDofValue('DISPLACEMENT_Y','end');
                expectedDisplacementY = A * cos(omega_D * time + theta) * exp(-delta * time);
                testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                    'Within', AbsoluteTolerance(1e-3)))
                
            end
        end
        
        function PODtest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %solve POD
            POD_model = SolverTests.create12DofModel;
            PODsolver = MORStrategy(POD_model,'POD','nPODmodes',8);
            coarse_sampling = linspace(1e-4,3.188,41);
            fine_sampling = linspace(0,1,120);
            PODsolver.initialize(coarse_sampling);
            PODsolver.solve(fine_sampling);
            disp = abs(POD_model.getNode(12).getDofValue('DISPLACEMENT_X','all'));
            load('tests/input_data/test_data.mat','POD_disp');
            
            testCase.assertThat(disp, IsEqualTo(POD_disp, ...
                    'Within', AbsoluteTolerance(1e-5)))
        end
        
        function blochInverse1DSolverTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            model = FemModel();
            
            %nodes
            model.addNewNode(1,0,0,0);
            for id = 1:3
                model.addNewNode(1+id,id,0,0);
            end
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', ...
                'ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
            
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('ROTATION_X');
            model.getAllNodes.fixDof('ROTATION_Z');
            
            %elements
            for id = 1:3
                model.addNewElement('BeamElement3d2n', id, [id id+1]);
            end
            model.getAllElements.setPropertyValue('IY',6.67e-4);
            model.getAllElements.setPropertyValue('IZ',6.67e-4);
            model.getAllElements.setPropertyValue('IT',6.67e-4);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
            model.getAllElements.setPropertyValue('CROSS_SECTION',0.2);
            model.getAllElements.setPropertyValue('DENSITY',2699);
            
            solver = BlochInverse1DSolvingStrategy(model);
            numberOfPhases = 60;
            numberOfBands = 4;
            
            [phases,frequencies]=solver.solve(numberOfPhases,numberOfBands);
            
            actual_f1 = frequencies(2,1); %first frequency of (quasi-)longitudinal wave
            actual_k1 = phases(1); %first phase
            l = 3; %length of the model
            rho = model.getElement(1).getPropertyValue('DENSITY');
            A = model.getElement(1).getPropertyValue('CROSS_SECTION');
            IY = model.getElement(1).getPropertyValue('IY');
            E = model.getElement(1).getPropertyValue('YOUNGS_MODULUS');
            
            actual_speedquasilongitudinalwave = (2*pi*l*actual_f1/actual_k1); %v = 2pi*l*delta_f/delta_k , second point is the origin (0/0)
            expected_speedquasilongitudinalwave = sqrt(E/rho);
            
            actual_f_at_pi = frequencies(1,numberOfPhases); %lowest frequency at phase pi
            expected_f = pi/(2*l^2)*sqrt(E*IY/(rho*A));
            
            testCase.assertThat(actual_speedquasilongitudinalwave, IsEqualTo(expected_speedquasilongitudinalwave, ...
                'Within', RelativeTolerance(1e-4)))
            
            testCase.assertThat(actual_f_at_pi, IsEqualTo(expected_f, ...
                'Within', RelativeTolerance(5e-3)))
        end
        
        function blochDirect1DSolverTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            io=MdpaInput('input_data/Federtests_Kreis.mdpa');
            model = io.readModel();
            model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
            
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e10);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
            model.getAllElements.setPropertyValue('DENSITY',2699);
            
            allNodes = model.getAllNodes();
            massNodeID = length(allNodes)+1;
            
            allElements = model.getAllElements();
            massElementID = length(allElements)+1;
            spring1ElementID = massElementID+1;
            spring2ElementID = massElementID+2;
            
            springNodes = model.getModelPart('GENERIC_0Grad').getNodes();
            
            leftSpringNode = springNodes(1,1);
            leftSNId = getId(leftSpringNode);
            
            rightSpringNode = springNodes(1,2);
            rightSNId = getId(rightSpringNode);
            
            model.addNewNode(massNodeID,0.1,0.075,0);
            model.addNewElement('SpringDamperElement3d2n',spring1ElementID,[leftSNId massNodeID]);
            model.addNewElement('SpringDamperElement3d2n',spring2ElementID,[massNodeID rightSNId]);
            model.addNewElement('ConcentratedMassElement3d1n',massElementID, massNodeID);
            
            model.getNode(leftSNId).addDof("DISPLACEMENT_Z");
            model.getNode(rightSNId).addDof("DISPLACEMENT_Z");
            model.getNode(massNodeID).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);
            
            model.getNode(leftSNId).fixDof('DISPLACEMENT_Z');
            model.getNode(rightSNId).fixDof('DISPLACEMENT_Z');
            model.getNode(massNodeID).fixDof('DISPLACEMENT_Z');
            model.getNode(massNodeID).fixDof('DISPLACEMENT_Y');
            
            model.getElement(spring1ElementID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
            model.getElement(spring1ElementID).setPropertyValue('ELEMENTAL_DAMPING',0);
            model.getElement(spring2ElementID).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
            model.getElement(spring2ElementID).setPropertyValue('ELEMENTAL_DAMPING',0);
            
            model.getElement(massElementID).setPropertyValue('ELEMENTAL_MASS',7e0);
            
            solver = BlochDirect1DSolvingStrategy(model);
            
            frequencies=linspace(0,5,10);
            [eigenValues]=solver.solve(frequencies);
            
            actualSolution = eigenValues(4,1:3);
            expectedSolution=[55208020.2355578 + 0.00000000000000i, 67661057.4948086 + 16658472.2286909i, 13310862.4339418 - 40288296.0224401i];
            
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', RelativeTolerance(1e-10)))
        end
        
        function prescribedDofSolverTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            [model, x0, ~, ~, ~] = createRectangularPlate(1, 1, 1, 1, 'elementType', 'QuadrilateralElement2d4n');
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
            
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',30e6);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0.0);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',3);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            x0.fixAllDofs();
            
            model.getNode(4).prescribeDof('DISPLACEMENT_Y',-1.55555555555556e-07);
            
            solver = HarmonicSolvingStrategy(model);
            solver.solve(0);
            
            actualDisplacement = model.getDofArray.getValue('end');
            expectedDisplacement=[0;0;-6.66666666666667e-08;...
                -1.11111111111111e-07;0;0;6.66666666666667e-08;-1.55555555555556e-07];
            
            testCase.assertThat(actualDisplacement, IsEqualTo(expectedDisplacement, ...
                'Within', RelativeTolerance(1e-10)))
        end
        
        function harmonicSolverTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            [model, x0, ~, ~, ~] = createRectangularPlate(1, .2, 2, 2, 'elementType', 'ShellElement3d4n');
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', ...
                'ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
            
            model.setPropertyValue('YOUNGS_MODULUS',2.1e11);
            model.setPropertyValue('POISSON_RATIO',0.3);
            model.setPropertyValue('DENSITY',7860);
            model.setPropertyValue('THICKNESS',0.001);
            model.setPropertyValue('NUMBER_GAUSS_POINT',2);
            model.addNewCondition('PointLoadCondition3d1n', 1, 5);
            model.getNode(5).setPropertyValue('POINT_LOAD',[0 0 -10]);
            x0.fixAllDofs();
            
            solver = HarmonicSolvingStrategy(model);
            freq = linspace(0.1,15,10);
            
            for ii=freq
                solver.solve(ii);
            end
            
            actualDisplacement = model.getNode(5).getDofValue('DISPLACEMENT_Z','all');
            expectedDisplacement = [0 -0.110398299437391 -0.124681359646738 ...
                -0.203788756268610 0.910555117679549 0.0815530371548063 ...
                0.0275266834778880 0.00806975947554387 -0.00238259577493668 ...
                -0.00957598196504360 -0.0156850242530438];
            
            testCase.assertThat(actualDisplacement, IsEqualTo(expectedDisplacement, ...
                'Within', RelativeTolerance(1e-10)))
        end
        
    end
    
    methods (Access = private, Static)
        function model = create12DofModel
            model = FemModel();
            
            %nodes
            n01 = model.addNewNode(1,0,0,0);
            n=12;
            for id = 1:n
                model.addNewNode(1+id,10*id/n,0,0);
            end
            endnode = model.getNode(n+1);
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            %elements
            for id = 1:n
                model.addNewElement('SpringDamperElement3d2n', id, [id id+1]);
            end
            model.getAllElements.setPropertyValue('ELEMENTAL_STIFFNESS',1);
            model.getAllElements.setPropertyValue('ELEMENTAL_DAMPING',0.6);
            
            %masses
            mass_id = length(model.getAllElements()) + 1;
            for id = 1:n
                model.addNewElement('ConcentratedMassElement3d1n',id + n, id+1);
            end
            
            model.getElements(mass_id:length(model.getAllElements)).setPropertyValue('ELEMENTAL_MASS',1);
            
            %boundary conditions
            addPointLoad(endnode,10,[1 0 0]);
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            n01.fixAllDofs();
        end
    end
    
end

