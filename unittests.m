classdef unittests < matlab.unittest.TestCase

    properties
        Model111
        Model243
        epsilon
    end
 
    methods(TestMethodSetup)
        function createaModel111(testCase)
            testCase.Model111 = createpde();
            W = 1;
            D = 1;
            H = 1;
            testCase.Model111.Geometry = multicuboid(W,D,H);
            generateMesh(testCase.Model111,'GeometricOrder','linear');
        end

        function createaModel243(testCase)
            testCase.Model243 = createpde();
            W = 2;
            D = 4;
            H = 3;
            testCase.Model243.Geometry = multicuboid(W,D,H);
            generateMesh(testCase.Model243,'GeometricOrder','linear');
        end

        function defineEpsilon(testCase)
            testCase.epsilon = 1E-6;
        end
    end
 
    methods(TestMethodTeardown)
    end
    
    methods (Test)
        function testWaterPlaneTriangulation1(testCase)
            [CL,~] = waterplaneTriangulation(testCase.Model111.Mesh,2);
            MTO = triangulation(testCase.Model111.Mesh.Elements.', testCase.Model111.Mesh.Nodes.');
            [MCL,~] = freeBoundary(MTO);
            [e,~] = size(CL);
            [Me,~] = size(MCL);

            testCase.verifyEqual(e*6,Me);
        end
        
        function testWaterPlaneTriangulation2(testCase)
            [CL,P] = waterplaneTriangulation(testCase.Model111.Mesh,2);
            MTO = triangulation(testCase.Model111.Mesh.Elements.', testCase.Model111.Mesh.Nodes.');
            [MCL,MP] = freeBoundary(MTO);

            testCase.verifyEqual(union(P(CL),MP(MCL)),unique(MP(MCL)));
        end

        function testWaterPlaneTriangulation3(testCase)
            [CL,P] = waterplaneTriangulation(testCase.Model111.Mesh,2);
            Z = P(unique(CL),3);
            [Ze,~] = size(Z);

            testCase.verifyLessThan(abs(sum(Z) - Ze),testCase.epsilon);
        end

        function testWaterPlaneTriangulation4(testCase)
            [CL,~] = waterplaneTriangulation(testCase.Model243.Mesh,2);
            [CL1,~] = waterplaneTriangulation(testCase.Model243.Mesh,1);
            [e,~] = size(CL);
            [e1,~] = size(CL1);

            testCase.verifyEqual(e,e1);
        end
        
        function testWaterPlaneTriangulation5(testCase)
            [CL,P] = waterplaneTriangulation(testCase.Model243.Mesh,2);
            MTO = triangulation(testCase.Model243.Mesh.Elements.', testCase.Model243.Mesh.Nodes.');
            [MCL,MP] = freeBoundary(MTO);

            testCase.verifyEqual(union(P(CL),MP(MCL)),unique(MP(MCL)));
        end

        function testWaterPlaneTriangulation6(testCase)
            [CL,P] = waterplaneTriangulation(testCase.Model243.Mesh,2);
            Z = P(CL,3);
            [Ze,~] = size(Z);

            testCase.verifyLessThan(abs(sum(Z) - 3*Ze),testCase.epsilon);
        end

        function testWaterPlaneTriangulationAndMoments1(testCase)
            [CL,P] = waterplaneTriangulation(testCase.Model111.Mesh,2);
            [S,~,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(S - 1),testCase.epsilon);
        end

        function testWaterPlaneTriangulationAndMoments2(testCase)
            [CL,P] = waterplaneTriangulation(testCase.Model111.Mesh,2);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [0;0;1]),testCase.epsilon);
        end

        function testWaterPlaneTriangulationAndMoments3(testCase)
            [CL,P] = waterplaneTriangulation(testCase.Model111.Mesh,3);
            [~,CF,~,~] = waterplaneMoments(CL,P);

            testCase.verifyLessThan(abs(CF - [0.5;0;0.5]),testCase.epsilon);
        end
    end
end

