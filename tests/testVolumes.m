classdef testVolumes < matlab.unittest.TestCase

    properties
        ModelSqr1
        ModelRect243
        ModelSphere1Off1
        ModelSphere2Offn2
        epsilon
    end    
 
    methods(TestMethodSetup)
        function createModelSqr1(testCase)
            testCase.ModelSqr1 = createpde();
            W = 1;
            D = 1;
            H = 1;
            testCase.ModelSqr1.Geometry = multicuboid(W,D,H);
            generateMesh(testCase.ModelSqr1,'GeometricOrder','linear');
        end

        function createModelRect243(testCase)
            testCase.ModelRect243 = createpde();
            W = 2;
            D = 4;
            H = 3;
            testCase.ModelRect243.Geometry = multicuboid(W,D,H);
            generateMesh(testCase.ModelRect243,'GeometricOrder','linear');
        end

        function createModelSphere1Off1(testCase)
            testCase.ModelSphere1Off1 = createpde();
            R = 1;
            Offset = [1,1,1];
            testCase.ModelSphere1Off1.Geometry = multisphere(R);
            translate(testCase.ModelSphere1Off1.Geometry, Offset);
            generateMesh(testCase.ModelSphere1Off1,'GeometricOrder','linear');
        end

        function createModelSphere2Offn2(testCase)
            testCase.ModelSphere2Offn2 = createpde();
            R = 2;
            Offset = -[2,2,2];
            testCase.ModelSphere2Offn2.Geometry = multisphere(R);
            translate(testCase.ModelSphere2Offn2.Geometry, Offset);
            generateMesh(testCase.ModelSphere2Offn2,'GeometricOrder','linear');
        end

        function defineEpsilon(testCase)
            testCase.epsilon = 1E-6;
        end
        
        function addPath(testCase)
            addpath("..\");
        end
    end
 
    methods(TestMethodTeardown)
        function removePath(testCase)
            rmpath("..\");
        end
    end
    
    methods(Test)
        function testVolumeMoments1(testCase)
            [V,~] = volumeMoments(testCase.ModelSqr1.Mesh);
            value = 1;

            testCase.verifyLessThan(abs(V-value),testCase.epsilon);
        end

        function testVolumeMoments2(testCase)
            [V,~] = volumeMoments(testCase.ModelRect243.Mesh);
            value = 24;

            testCase.verifyLessThan(abs(V-value),testCase.epsilon);
        end
        
        function testVolumeMoments3(testCase)
            [V,~] = volumeMoments(testCase.ModelSphere1Off1.Mesh);
            value = 4/3*pi*(1)^3;

            testCase.verifyLessThan(abs(V-value),0.1);
        end

        function testVolumeMoments4(testCase)
            [V,~] = volumeMoments(testCase.ModelSphere2Offn2.Mesh);
            value = 4/3*pi*(2)^3;

            testCase.verifyLessThan(abs(V-value),0.5);
        end

         function testVolumeMoments5(testCase)
            [~,CB] = volumeMoments(testCase.ModelSqr1.Mesh);
            value = [0;0;0.5];

            testCase.verifyLessThan(abs(CB-value),testCase.epsilon);
        end

        function testVolumeMoments6(testCase)
            [~,CB] = volumeMoments(testCase.ModelRect243.Mesh);
            value = [0;0;1.5];

            testCase.verifyLessThan(abs(CB-value),testCase.epsilon);
        end
        
        function testVolumeMoments7(testCase)
            [~,CB] = volumeMoments(testCase.ModelSphere1Off1.Mesh);
            value = [1;1;1];

            testCase.verifyLessThan(abs(CB-value),5E-3);
        end

        function testVolumeMoments8(testCase)
            [~,CB] = volumeMoments(testCase.ModelSphere2Offn2.Mesh);
            value = -[2;2;2];

            testCase.verifyLessThan(abs(CB-value),5E-3);
        end
    end
end