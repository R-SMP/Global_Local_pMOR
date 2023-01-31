classdef ValidationTests <  matlab.unittest.TestCase
    %VALIDATIONTESTS Class for larger tests
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        % bridge test taken from IFEM ch. 21 by Felippa
        function bridgeTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            io = GmshInput('validation_bridge_input.msh');
            model = io.readModel;
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            model.getModelPart('fixed_support').getNodes.fixAllDofs;
            model.getModelPart('roller_support').getNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            
            addPointLoad(model.getNodes([3 5 9 11]),10,[0 -1 0]);
            addPointLoad(model.getNode(7),16,[0 -1 0]);
            
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacementX = model.getAllNodes.getDofValue('DISPLACEMENT_X');
            actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
            
            expectedDisplacementX = [0 0.809536 0.28 0.899001 0.56 0.8475 ...
                0.8475 0.795999 1.135 0.885464 1.415 1.695]';
            expectedDisplacementY = [0 -1.775600 -1.792260 -2.291930 -2.316600 ...
                -2.385940 -2.421940 -2.291930 -2.316600 -1.775600 -1.792260 0]';
            
            testCase.assertThat(actualDisplacementX, IsEqualTo(expectedDisplacementX, ...
                'Within', RelativeTolerance(1e-5)))
            testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                'Within', RelativeTolerance(1e-5)))
            
            actualElementStress = model.getAllElements.computeElementStress(1);
            expectedElementStress = [28 28 28.75 28.75 28 28 -6.2610 -6.0030 ...
                -6.0300 -6.0300 -6.0030 -6.2610 3.3330 3.0830 4.0000 3.0830 ...
                3.3330 1.6770 3.2020 3.2020 1.6770]';
            
            testCase.assertThat(actualElementStress, IsEqualTo(expectedElementStress, ...
                'Within', RelativeTolerance(1e-3)))
        end
        
        function sdofWithHarmonicExcitation(testCase)
            %SDOFWITHHARMONICEXCITATION sdof spring mass system with
            %harmonic excitiation
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %parameters
            mass = 80.0;
            stiffness = 200.0;
            F0 = 100.0;
            uInit = 0.0;
            vInit = 0.0;
            
            %analytical solution
            eigenfrequency = sqrt(stiffness/mass);
            excitationFrequency = 0.5 * eigenfrequency;
            beta = excitationFrequency / eigenfrequency;
            
            model = FemModel();
            n01 = model.addNewNode(1,0,0,0);
            n02 = model.addNewNode(2,0,1,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            s = model.addNewElement('SpringDamperElement3d2n',1,[1 2]);
            s.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            m = model.addNewElement('ConcentratedMassElement3d1n',2,2);
            m.setPropertyValue('ELEMENTAL_MASS',mass);
            
            model.getAllNodes.fixDof('DISPLACEMENT_X');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            n01.fixDof('DISPLACEMENT_Y');
            
            dt = 0.05;
            time = 0;
            endTime = 20;
            solver = NewmarkSolvingStrategy(model, dt);
            while time < endTime
                applyHarmonicSineLoad(n02, F0, [0, -1, 0], excitationFrequency, time);
                solver.solve();
                
                actualDisplacementY = n02.getDofValue('DISPLACEMENT_Y','end');
                transientY = uInit * cos(eigenfrequency * time) + vInit / eigenfrequency * ...
                    sin(eigenfrequency * time) - F0 / stiffness * beta / (1 - power(beta,2)) * ...
                    sin(eigenfrequency * time);
                steadyY = F0 / stiffness / (1 - power(beta, 2)) * sin(excitationFrequency * time);
                expectedDisplacementY = - (transientY + steadyY);
                testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                    'Within', AbsoluteTolerance(1e-2)))
                
                time = time + dt;
            end
            
        end
        
        function mdofWithDamping(testCase)
            %MDOFWITHDAMPING 2-dof-system with springs, masses, and dampers;
            %taken from HUMAR: Dynamics for Structures (p. 560)
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.RelativeTolerance
            tolerance = AbsoluteTolerance(0.07) | RelativeTolerance(0.07);
            
            model = FemModel();
            n01 = model.addNewNode(1,1,0,0);
            n02 = model.addNewNode(2,2,0,0);
            base = model.addNewNode(3,0,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            s01 = model.addNewElement('SpringDamperElement3d2n',1,[3 1]);
            s01.setPropertyValue('ELEMENTAL_STIFFNESS', 2.0);
            s01.setPropertyValue('ELEMENTAL_DAMPING', 0.4);
            m01 = model.addNewElement('ConcentratedMassElement3d1n',2,1);
            m01.setPropertyValue('ELEMENTAL_MASS', 2.0);
            s02 = model.addNewElement('SpringDamperElement3d2n',3,[1 2]);
            s02.setPropertyValue('ELEMENTAL_STIFFNESS', 1.0);
            s02.setPropertyValue('ELEMENTAL_DAMPING', 0.05);
            m02 = model.addNewElement('ConcentratedMassElement3d1n',4,2);
            m02.setPropertyValue('ELEMENTAL_MASS', 1.0);
            s03 = model.addNewElement('SpringDamperElement3d2n',5,[3 2]);
            s03.setPropertyValue('ELEMENTAL_DAMPING', 0.15);
            
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            base.fixDof('DISPLACEMENT_X');
            
            n01.setDofValue('DISPLACEMENT_X', 1.0);
            n02.setDofValue('DISPLACEMENT_X', 2.0);
            
            dt = 0.05;
            time = 0;
            endTime = 20;
            solver = NewmarkSolvingStrategy(model, dt);
            while time < endTime
                solver.solve();
                
                actualDisplacementX01 = n01.getDofValue('DISPLACEMENT_X','end');
                actualDisplacementX02 = n02.getDofValue('DISPLACEMENT_X','end');
                
                expectedDisplacementX01 = 1.0064 * exp(-0.08334 * time) * sin(0.70221 * time + atan(9.031));
                expectedDisplacementX02 = 2.0148 * exp(-0.08334 * time) * sin(0.70221 * time + atan(8.153));
                
                testCase.assertThat(actualDisplacementX01, IsEqualTo(expectedDisplacementX01, ...
                    'Within', tolerance))
                testCase.assertThat(actualDisplacementX02, IsEqualTo(expectedDisplacementX02, ...
                    'Within', tolerance))
                
                time = time + dt;
            end
            
        end
        
        function twoDofHarmonicAnalysis(testCase)
            %TWODOFHARMONICANALYSIS 2-dof-system with springs and masses;
            %taken from HUMAR: Dynamics for Structures (p. 675)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            model = FemModel();
            support = model.addNewNode(3,0,0,0);
            n01 = model.addNewNode(1,10,0,0);
            model.addNewNode(2,20,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            stiffness = 20.0;
            mass = 0.2;
            
            s01 = model.addNewElement('SpringDamperElement3d2n',1,[1 3]);
            s01.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            s02 = model.addNewElement('SpringDamperElement3d2n',2,[1 2]);
            s02.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness/2);
            m01 = model.addNewElement('ConcentratedMassElement3d1n',3,1);
            m01.setPropertyValue('ELEMENTAL_MASS', mass);
            m02 = model.addNewElement('ConcentratedMassElement3d1n',4,2);
            m02.setPropertyValue('ELEMENTAL_MASS', mass/2);
            
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            support.fixDof('DISPLACEMENT_X');
            
            addPointLoad(n01,1,[1 0 0]);
            
            exfreq = linspace(.1*sqrt(5),10*sqrt(5),1000);
            
            solver = EigensolverStrategy(model);
            solver.harmonicAnalysis(exfreq,2);
            
            Rd11_actual = model.getNode(1).getDofValue('DISPLACEMENT_X','all');
            Rd12_actual = model.getNode(2).getDofValue('DISPLACEMENT_X','all');
            
            Rd11_expected = (1/3 * stiffness) ./ (0.5 - (exfreq ./ sqrt(stiffness/mass)).^2) ...
                + (1/1.5 * stiffness) ./ (2 - (exfreq ./ sqrt(stiffness/mass)).^2);
            Rd12_expected = (2/3 * stiffness) ./ (0.5 - (exfreq ./ sqrt(stiffness/mass)).^2) ...
                - (2/3 * stiffness) ./ (2 - (exfreq ./ sqrt(stiffness/mass)).^2);
            
            testCase.assertThat(Rd11_actual, IsEqualTo(Rd11_expected / stiffness^2, ...
                'Within', AbsoluteTolerance(1e-7)))
            testCase.assertThat(Rd12_actual, IsEqualTo(Rd12_expected / stiffness^2, ...
                'Within', AbsoluteTolerance(1e-7)))
            
        end
        
        function shellElement3d4nLargeTest(testCase)
            %SHELLELEMENT3D4NLARGETEST plane shell consisting out of 100
            %   elements under static loading at the middle node
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 10, 10, 'elementType', 'ShellElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
            
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
            
            support = [x0 xl y0 yl];
            support.fixAllDofs();
            
            model.getNode(61).setDofLoad('DISPLACEMENT_Z',2500);
            
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacement = model.getNode(61).getDofValue('DISPLACEMENT_Z');
            expectedDisplacement = 0.006040637455055775;
            
            testCase.assertThat(actualDisplacement, IsEqualTo(expectedDisplacement, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function planeStressElement2d4nLargeTest(testCase)
            %PLANESTRESSELEMENT3D4NLARGETEST 160 4-node plane stress
            %   elements cantilevered with single load
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            io = GmshInput('tests/input_data/Quad4_4x40.msh');
            model = io.readModel;
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
            
            model.getAllElements.setPropertyValue('THICKNESS', 0.5);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT', 3);
            model.getAllElements.setPropertyValue('DENSITY', 7850);
            
            idLowerLCorner = model.getModelPart('LLCorner').getNodes().getId();
            idUpperLCorner = model.getModelPart('ULCorner').getNodes().getId();
            addLineBC(idLowerLCorner,idUpperLCorner,model.getAllNodes);
            
            model.getModelPart('URCorner').getNodes().setDofLoad('DISPLACEMENT_Y',-1000);
            
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacement = model.getModelPart('URCorner').getNodes().getDofValue('DISPLACEMENT_Y');
            expectedDisplacement = -3.90385982503743e-05;
            
            testCase.assertThat(actualDisplacement, IsEqualTo(expectedDisplacement, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function planeStressElement2d6nLargeTest(testCase)
            %PLANESTRESSELEMENT3D6NLARGETEST 544 6-node plane stress
            %   elements cantilevered with single load
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            io = MdpaInput('tests/input_data/ecke01.mdpa');
            m = io.readModel();
            
            m.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
            
            m.getAllElements.setPropertyValue('THICKNESS', 0.5);
            m.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2e11);
            m.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            m.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT', 3);
            m.getAllElements.setPropertyValue('DENSITY', 7850);
            
            m.getModelPart('GENERIC_support').getNodes().fixAllDofs();
            m.getModelPart('GENERIC_load').getNodes().setDofLoad('DISPLACEMENT_Y',-1e7);
            
            solver = SimpleSolvingStrategy(m);
            solver.solve();
            
            actualDisplacement = m.getNode(81).getDofValue('DISPLACEMENT_Y');
            expectedDisplacement = -0.0496445026676462;
            
            testCase.assertThat(actualDisplacement, IsEqualTo(expectedDisplacement, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function shellElementTransformationTest(testCase)
            %SHELLELEMENTTRANSFORMATIONTEST A plate of 100 ShellElement3d4n
            %   not parallel to the xyz coordinate system
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            % read model
            io = MdpaInput('tests/input_data/shell_transformation_test.mdpa');
            m = io.readModel();
            m.getAllNodes.addDof(["DISPLACEMENT_X", "DISPLACEMENT_Y", ...
                "DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
            
            % get normal vector
            e = m.getElement(1);
            en = e.getNodes();
            n = cross(en(1).getCoords(),en(4).getCoords());
            n = n/norm(n);
            load = 1000*n;
            
            % material properties
            m.setPropertyValue('THICKNESS', 0.005);
            m.setPropertyValue('YOUNGS_MODULUS', 206.9e9);
            m.setPropertyValue('POISSON_RATIO', 0.29);
            m.setPropertyValue('NUMBER_GAUSS_POINT', 2);
            m.setPropertyValue('DENSITY', 7850);
            
            % boundary conditions and loads
            m.getModelPart('GENERIC_support').getNodes().fixAllDofs();
            m.getModelPart('GENERIC_load').getNodes().setPropertyValue('POINT_LOAD',load);
            
            % solution
            solver = SimpleSolvingStrategy(m);
            solver.solve();
            
            % assertion
            actual_x = m.getModelPart('GENERIC_load').getNodes().getDofValue('DISPLACEMENT_X');
            actual_y = m.getModelPart('GENERIC_load').getNodes().getDofValue('DISPLACEMENT_Y');
            actual_z = m.getModelPart('GENERIC_load').getNodes().getDofValue('DISPLACEMENT_Z');
            actual_abs = sqrt(actual_x^2+actual_y^2+actual_z^2);
            
            testCase.assertThat(actual_x, IsEqualTo(0.00223774381522016, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actual_y, IsEqualTo(0.000241366892548215, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actual_z, IsEqualTo(-0.001006380598812, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actual_abs, IsEqualTo(0.002465473031497, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function BlochInverse1DValidationTest(testCase)
            % function Validation_Test_BlochInverse1D()
            %%% unit cell with l = 0.2m, h = 0.15m; a round inclusion in the center
            %%% with r = 0.025m; a horizontal spring mass system will be added in the
            %%% inclusion
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            io=MdpaInput('tests/input_data/Federtests_Kreis.mdpa'); %specify input file
            model = io.readModel(); %read the model
            
            model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
            
            %Element properties for homogenous unic-cell: here aluminium
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
            
            % Create Bloch Solver
            solver = BlochInverse1DSolvingStrategy(model);
            % define number of phases and number of bands
            numberOfPhases = 50;
            numberOfBands = 10;
            
            [phases,frequencies]=solver.solve(numberOfPhases,numberOfBands);
            
            actual_f_longitudinal = frequencies(2,1);  %first frequency of (quasi-)longitudinal wave
            actual_phase_longitudinal = phases(1);  %checks starting point of phases
            actual_f_bending_at_pi = frequencies(2,numberOfPhases); %Frequency of bending wave at phase pi
            
            expected_f_longitudinal = 36.469358906958256;
            expected_f_bending_at_pi = 4.502029917504889e+03;
            expected_phase_longitudinal = 0.01;
            
            testCase.assertThat(actual_f_longitudinal, IsEqualTo(expected_f_longitudinal, ...
                'Within', RelativeTolerance(1e-4)))
            
            testCase.assertThat(actual_f_bending_at_pi, IsEqualTo(expected_f_bending_at_pi, ...
                'Within', RelativeTolerance(1e-4)))
            
            testCase.assertThat(actual_phase_longitudinal, IsEqualTo(expected_phase_longitudinal, ...
                'Within', RelativeTolerance(1e-4)))  %if this doesn't work, in the
            %BlochInverse1D solver the static method 'propConst(nOP)' was
            %changed
            
        end
        
        function externalScriptsTest(testCase)
            bridge;
            Bridge_with_inputFile;
        end
        
        function testSoilSuperElement(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            Nx = 64;
            dx = 0.3;

            nodeArray = generateSoilSuperElement(dx,Nx);
            
            nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

            ele01 = SoilSuperElement(1,nodeArray,Nx);

            % properties of upper layer
            ele01.setPropertyValue('THICKNESS',1);                              
            ele01.setPropertyValue('YOUNGS_MODULUS',45e6);
            ele01.setPropertyValue('POISSON_RATIO',0.332);
            ele01.setPropertyValue('DENSITY',1800);
            ele01.setPropertyValue('ELEMENTAL_DAMPING',0.05);

            % properties of half space
            ele01.setPropertyValue('YOUNGS_MODULUS_HS',45e6);        
            ele01.setPropertyValue('POISSON_RATIO_HS',0.332);
            ele01.setPropertyValue('DENSITY_HS',1800);
            ele01.setPropertyValue('ELEMENTAL_DAMPING_HS',0.05);
            
            % discretization parameters
            ele01.setPropertyValue('SAMPLE_NUMBER',Nx);          
            ele01.setPropertyValue('SAMPLE_RATE',dx);
            ele01.setPropertyValue('FREQUENCY',25);
            
            % point load 
            j = getLocalNodeId(nodeArray,0,0);
            addPointLoad(nodeArray(j),6,[1 0 0]);
            
            % model
            model = FemModel(nodeArray, ele01);
            
            % solver
            solver = SimpleSolvingStrategy(model);
            solver.solve(); 
            
            % check results at node (x,y) = (dx,dx)
            j = getLocalNodeId(nodeArray,dx,dx);
            check_node = model.getNode(j);
            actualSolution = [check_node.getDofValue('DISPLACEMENT_X',1);...
                              check_node.getDofValue('DISPLACEMENT_Y',1);...
                              check_node.getDofValue('DISPLACEMENT_Z',1)];
            
            expectedSolution = [6.6936e-07 + 5.6376e-07i;
                                2.6542e-07 + 3.2420e-08i;
                                1.6486e-07 + 7.7518e-08i];
            
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within',AbsoluteTolerance(1e-11)))
        end
    end
    
end
