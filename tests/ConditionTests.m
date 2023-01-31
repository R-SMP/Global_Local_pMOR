classdef ConditionTests < matlab.unittest.TestCase
    %ELEMENTTESTS Tests for all elements
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        function testAcousticNeumannLineCondition2d2n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            io = MdpaInput('Acoustic_domain_fixed.mdpa');
            model = io.readModel;
            
            model.getAllNodes.addDof({'ACOUSTIC_PRESSURE'});
            
            neumann_nodes = model.getModelPart('LinePressure2D_Boundary').getNodes.getId;
            neumann_condition = model.getModelPart('LinePressure2D_Boundary').getConditions.getId;
            
            for i=1:length(neumann_condition)
                model.getCondition(neumann_condition(i)).setPropertyValue('FREQUENCY',200*2*pi);
            end
            for i=1:length(neumann_nodes)
                model.getNode(neumann_nodes(i)).setPropertyValue('NORMAL_DISPLACEMENT',0);
            end
            
            model.getNode(15).prescribeDof('ACOUSTIC_PRESSURE',1);
            

            solver = HarmonicSolvingStrategy(model);
            solver.solve(200*2*pi);
            
            actualSolution = model.getNode(25).getDof('ACOUSTIC_PRESSURE').getValue('end');
            
            expectedSolution = -2.52717993554755;
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', AbsoluteTolerance(1e-7)))
        end
        
        function testAcousticRobinLineCondition2d2n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            io = MdpaInput('Acoustic_domain_free.mdpa');
            model = io.readModel;
            
            model.getAllNodes.addDof({'ACOUSTIC_PRESSURE'});
            
            robin_nodes = model.getModelPart('LinePressure2D_Boundary').getNodes.getId;
            
            for i=1:length(robin_nodes)
                model.getNode(robin_nodes(i)).setPropertyValue('ADMITTANCE',2.41e-3);
            end
            
            model.getNode(15).prescribeDof('ACOUSTIC_PRESSURE',1);
            
            solver = HarmonicSolvingStrategy(model);
            solver.solve(10*2*pi);
            
            actualSolution = model.getNode(25).getDof('ACOUSTIC_PRESSURE').getValue('end');
            
            expectedSolution = 0.897407226090890 - 0.313206884845964i;
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', AbsoluteTolerance(1e-7)))
        end
        
        function testAcousticStructureLineCondition2d2n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            io = MdpaInput('Fluid_Structure.mdpa');
            model = io.readModel;
            
            Fluid = model.getModelPart('Parts_Acoustic').getNodes.getId;
            for i=1:length(Fluid)
                model.getNode(Fluid(i)).addDof({'ACOUSTIC_PRESSURE'});
            end
            
            Structure = model.getModelPart('Parts_Structure').getNodes.getId;
            for i=1:length(Structure)
                model.getNode(Structure(i)).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
            end
            
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',2.1e6);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
            model.getAllElements.setPropertyValue('DENSITY',250);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
            
            neumann_nodes = model.getModelPart('LinePressure2D_Neumann').getNodes.getId;
            neumann_condition = model.getModelPart('LinePressure2D_Neumann').getConditions.getId;
            for i=1:length(neumann_condition)
                model.getCondition(neumann_condition(i)).setPropertyValue('FREQUENCY',200*2*pi);
            end
            for i=1:length(neumann_nodes)
                model.getNode(neumann_nodes(i)).setPropertyValue('NORMAL_DISPLACEMENT',0);
            end
            
            Side = model.getModelPart('GENERIC_Side').getNodes.getId;
            for i=1:length(Side)
                model.getNode(Side(i)).prescribeDof('DISPLACEMENT_X',0);
                model.getNode(Side(i)).prescribeDof('DISPLACEMENT_Y',0);
            end
            
            model.getNode(4).prescribeDof('ACOUSTIC_PRESSURE',1);
            
            solver = HarmonicSolvingStrategy(model);
            solver.solve(200*2*pi);
            
            actualSolution = model.getNode(11).getDof('DISPLACEMENT_X').getValue('end');
            
            expectedSolution = 1.09484735787404e-08;
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', AbsoluteTolerance(1e-7)))
        end
        
        function testAcousticCPorousLineCondition2d2n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            io = MdpaInput('ClassicalPorous_Acoustic.mdpa');
            model = io.readModel;
            
            Fluid = model.getModelPart('Parts_Acoustic').getNodes.getId;
            for i=1:length(Fluid)
                model.getNode(Fluid(i)).addDof({'ACOUSTIC_PRESSURE'});
            end
            
            Porous = model.getModelPart('Parts_Porous').getNodes.getId;
            for i=1:length(Porous)
                model.getNode(Porous(i)).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_FLUID_X","DISPLACEMENT_FLUID_Y"]);
            end
            
            model.getAllElements.setPropertyValue('DENSITY_SOLID',750);
            model.getAllElements.setPropertyValue('LAMBDA_SOLID',905357);
            model.getAllElements.setPropertyValue('MUE_SOLID',264062);
            model.getAllElements.setPropertyValue('DAMPING_SOLID',0);
            model.getAllElements.setPropertyValue('POROSITY',0.96);
            model.getAllElements.setPropertyValue('TORTUOSITY',1.7);
            model.getAllElements.setPropertyValue('FLOW_RESISTIVITY',32e3);
            model.getAllElements.setPropertyValue('VISCOUS_LENGHT',226e-6);
            model.getAllElements.setPropertyValue('THERMAL_LENGTH',226e-6);
            model.getAllElements.setPropertyValue('THERMAL_LENGTH',226e-6);
            model.getAllElements.setPropertyValue('FREQUENCY',2*pi*100);
            
            neumann_nodes = model.getModelPart('LinePressure2D_Neumann').getNodes.getId;
            neumann_condition = model.getModelPart('LinePressure2D_Neumann').getConditions.getId;
            for i=1:length(neumann_condition)
                model.getCondition(neumann_condition(i)).setPropertyValue('FREQUENCY',100*2*pi);
            end
            for i=1:length(neumann_nodes)
                model.getNode(neumann_nodes(i)).setPropertyValue('NORMAL_DISPLACEMENT',0);
            end
            
            interface_condition = model.getModelPart('LinePressure2D_Interface').getConditions.getId;
            for i=1:length(interface_condition)
                model.getCondition(interface_condition(i)).setPropertyValue('POROSITY',0.96);
            end
            
            Side = model.getModelPart('GENERIC_Side').getNodes.getId;
            for i=1:length(Side)
                model.getNode(Side(i)).prescribeDof('DISPLACEMENT_Y',0);
                model.getNode(Side(i)).prescribeDof('DISPLACEMENT_FLUID_Y',0);
            end
            
            Back = model.getModelPart('GENERIC_Back').getNodes.getId;
            for i=1:length(Back)
                model.getNode(Back(i)).prescribeDof('DISPLACEMENT_Y',0);
                model.getNode(Back(i)).prescribeDof('DISPLACEMENT_X',0);
                model.getNode(Back(i)).prescribeDof('DISPLACEMENT_FLUID_Y',0);
                model.getNode(Back(i)).prescribeDof('DISPLACEMENT_FLUID_X',0);
            end
            
            model.getNode(4).prescribeDof('ACOUSTIC_PRESSURE',1);
            
            solver = HarmonicSolvingStrategy(model);
            solver.solve(100*2*pi);
            
            actualSolution = model.getNode(11).getDof('DISPLACEMENT_X').getValue('end');
            
            expectedSolution = 1.15024509096351e-06 - 2.55260073901526e-07i;
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', AbsoluteTolerance(1e-7)))
        end
        
        function testStructureMPorousLineCondition2d2n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            io = MdpaInput('MixedPorous_Structure.mdpa');
            model = io.readModel;
            
            omega = 2*pi*100;
            
            Solid = model.getModelPart('Parts_Solid').getNodes.getId;
            for i=1:length(Solid)
                model.getNode(Solid(i)).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
            end
            
            Porous = model.getModelPart('Parts_Porous').getNodes.getId;
            for i=1:length(Porous)
                model.getNode(Porous(i)).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","ACOUSTIC_PRESSURE"]);
            end
            
            model.getAllElements.setPropertyValue('DENSITY_SOLID',750);
            model.getAllElements.setPropertyValue('LAMBDA_SOLID',905357);
            model.getAllElements.setPropertyValue('MUE_SOLID',264062);
            model.getAllElements.setPropertyValue('DAMPING_SOLID',0);
            model.getAllElements.setPropertyValue('POROSITY',0.96);
            model.getAllElements.setPropertyValue('TORTUOSITY',1.7);
            model.getAllElements.setPropertyValue('FLOW_RESISTIVITY',32e3);
            model.getAllElements.setPropertyValue('VISCOUS_LENGHT',90e-6);
            model.getAllElements.setPropertyValue('THERMAL_LENGTH',165e-6);
            
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',2.1e6);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
            model.getAllElements.setPropertyValue('DENSITY',250);
            
            model.getAllElements.setPropertyValue('FREQUENCY',omega);
            
            Back = model.getModelPart('GENERIC_Back').getNodes.getId;
            Front = model.getModelPart('GENERIC_Front').getNodes.getId;
            Side = model.getModelPart('GENERIC_Side').getNodes.getId;
            for i=1:length(Side)
                model.getNode(Side(i)).prescribeDof('DISPLACEMENT_X',0);
            end
            for i=1:length(Back)
                model.getNode(Back(i)).prescribeDof('DISPLACEMENT_X',0);
                model.getNode(Back(i)).prescribeDof('DISPLACEMENT_Y',0);
            end
            for i=1:length(Front)
                model.getNode(Front(i)).prescribeDof('DISPLACEMENT_Y',1e-6);
            end
            
            solver = HarmonicSolvingStrategy(model);
            solver.solve(omega);
            
            actualSolution = model.getNode(180).getDof('ACOUSTIC_PRESSURE').getValue('end');
            
            expectedSolution = -0.383898658433232 - 0.0304939371188891i;
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', AbsoluteTolerance(1e-7)))
        end
        
        function testSeptumMPorousLineCondition2d2n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            io = MdpaInput('MixedPorous_Septum.mdpa');
            model = io.readModel;
            
            omega = 2*pi*100;
            
            Porous = model.getModelPart('Parts_Porous').getNodes.getId;
            for i=1:length(Porous)
                model.getNode(Porous(i)).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","ACOUSTIC_PRESSURE"]);
            end
            
            model.getAllElements.setPropertyValue('DENSITY_SOLID',750);
            model.getAllElements.setPropertyValue('LAMBDA_SOLID',487500);
            model.getAllElements.setPropertyValue('MUE_SOLID',325000);
            model.getAllElements.setPropertyValue('DAMPING_SOLID',0.0);
            
            model.getAllElements.setPropertyValue('POROSITY',0.96);
            model.getAllElements.setPropertyValue('TORTUOSITY',1.7);
            model.getAllElements.setPropertyValue('FLOW_RESISTIVITY',32e3);
            model.getAllElements.setPropertyValue('VISCOUS_LENGHT',90e-6);
            model.getAllElements.setPropertyValue('THERMAL_LENGTH',165e-6);
            
            model.getAllElements.setPropertyValue('FREQUENCY',omega);
            
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',30e6);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',3);
            model.getAllElements.setPropertyValue('DENSITY',860);
            
            Support = model.getModelPart('GENERIC_Support').getNodes.getId;
            Excitation = model.getModelPart('GENERIC_Excitation').getNodes.getId;
            septum_condition = model.getModelPart('LinePressure2D_Septum').getConditions.getId;
            for i =1:length(Support)
                model.getNode(Support(i)).fixDof('DISPLACEMENT_X');
                model.getNode(Support(i)).fixDof('DISPLACEMENT_Y');
            end
            for i =1:length(Excitation)
                model.getNode(Excitation(i)).prescribeDof('ACOUSTIC_PRESSURE',1);
            end
            for i=1:length(septum_condition)
                model.getCondition(septum_condition(i)).setPropertyValue('SURFACE_DENSITY',6);
            end
            
            solver = HarmonicSolvingStrategy(model);
            solver.solve(omega);
            
            actualSolution = model.getNode(44).getDof('ACOUSTIC_PRESSURE').getValue('end');
            
            expectedSolution = 1.00011132969214 - 0.00872872195778083i;
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', AbsoluteTolerance(1e-7)))
        end
        
    end
end
