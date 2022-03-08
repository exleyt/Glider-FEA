classdef testWaterplanes < matlab.unittest.TestCase

    properties
        ModelSqr1
        ModelRect243
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
        function testWaterplaneTriangulation1(testCase)
            [CL,~] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,2);
            MTO = triangulation(testCase.ModelSqr1.Mesh.Elements.', testCase.ModelSqr1.Mesh.Nodes.');
            [MCL,~] = freeBoundary(MTO);
            [e,~] = size(CL);
            [Me,~] = size(MCL);

            testCase.verifyEqual(e*6,Me);
        end
        
        function testWaterplaneTriangulation2(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,2);
            MTO = triangulation(testCase.ModelSqr1.Mesh.Elements.', testCase.ModelSqr1.Mesh.Nodes.');
            [MCL,MP] = freeBoundary(MTO);

            testCase.verifyEqual(union(P(CL),MP(MCL)),unique(MP(MCL)));
        end

        function testWaterplaneTriangulation3(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,2);
            Z = P(unique(CL),3);
            [Ze,~] = size(Z);

            testCase.verifyLessThan(abs(sum(Z) - Ze),testCase.epsilon);
        end

        function testWaterplaneTriangulation4(testCase)
            [CL,~] = waterplaneTriangulation(testCase.ModelRect243.Mesh,2);
            [CL1,~] = waterplaneTriangulation(testCase.ModelRect243.Mesh,1);
            [e,~] = size(CL);
            [e1,~] = size(CL1);

            testCase.verifyEqual(e,e1);
        end
        
        function testWaterplaneTriangulation5(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,2);
            MTO = triangulation(testCase.ModelRect243.Mesh.Elements.', testCase.ModelRect243.Mesh.Nodes.');
            [MCL,MP] = freeBoundary(MTO);

            testCase.verifyEqual(union(P(CL),MP(MCL)),unique(MP(MCL)));
        end

        function testWaterplaneTriangulation6(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,2);
            Z = P(CL,3);
            [Ze,~] = size(Z);

            testCase.verifyLessThan(abs(sum(Z) - 3*Ze),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments1(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,2);
            [S,~,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(S - 1),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments2(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,2);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [0;0;0]),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments3(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,3);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [0.5;0;0]),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments4(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,6);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [0;-0.50;0]),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments5(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,2);
            [~,~,Ixx,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(Ixx-1/12),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments6(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelSqr1.Mesh,2);
            [~,~,~,Iyy] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(Iyy-1/12),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments7(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,2);
            [S,~,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(S - 8),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments8(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,3);
            [S,~,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(S - 12),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments9(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,4);
            [S,~,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(S - 6),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments10(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,2);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [0;0;0]),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments11(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,4);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [0;2;0]),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments12(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,5);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [-1;0;0]),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments13(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,2);
            [~,~,Ixx,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(Ixx-8/3),testCase.epsilon);
        end

        function testWaterplaneTriangulationAndMoments14(testCase)
            [CL,P] = waterplaneTriangulation(testCase.ModelRect243.Mesh,2);
            [~,~,~,Iyy] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(Iyy-32/3),testCase.epsilon);
        end
    end
end

