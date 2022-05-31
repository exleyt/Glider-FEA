classdef testWaterplanes < matlab.unittest.TestCase

    properties
        MeshCube
        ModelSquare
        ModelSquare2
        epsilon
    end
 
    methods(TestMethodSetup)
        function createModelCube(testCase)
            modelCube = createpde();
            W = 1;
            D = 1;
            H = 1;
            modelCube.Geometry = multicuboid(W,D,H);
            generateMesh(modelCube,'GeometricOrder','linear');
            TO = triangulation(modelCube.Mesh.Elements.', modelCube.Mesh.Nodes.');
            [CL,P] = freeBoundary(TO);
            testCase.MeshCube = triangulation(CL,P);
        end

        function createModelSquare(testCase)
            P = [0,0,0;         % 1
                 1,0,0;         % 2
                 1,1,0;         % 3
                 0,1,0;         % 4
                 0.5,0.5,0];    % 5
            CL = [1,2,5;
                  2,3,5;
                  3,4,5;
                  4,1,5];
            testCase.ModelSquare = triangulation(CL,P);
        end

        function createModelSquare2(testCase)
            P = [1,1,-3;        % 1
                 4,2,-3;        % 2
                 4,5,-3;        % 3
                 1,4,-3;        % 4
                 2,2,-3];       % 5
            CL = [1,2,5;
                  2,3,5;
                  3,4,5;
                  4,1,5];
            testCase.ModelSquare2 = triangulation(CL,P);
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
        function testWaterplaneMoments1(testCase)
            FN = -faceNormal(testCase.MeshCube);
            [S,~,~,~,~,~] = waterplaneMoments(testCase.MeshCube.ConnectivityList,testCase.MeshCube.Points,FN);

            testCase.verifyLessThan(abs(S),testCase.epsilon);
        end

        function testWaterplaneMoments2(testCase)
            FN = -faceNormal(testCase.MeshCube);
            [~,Sx,~,~,~,~] = waterplaneMoments(testCase.MeshCube.ConnectivityList,testCase.MeshCube.Points,FN);

            testCase.verifyLessThan(abs(Sx),testCase.epsilon);
        end

        function testWaterplaneMoments3(testCase)
            FN = -faceNormal(testCase.MeshCube);
            [~,~,Sy,~,~,~] = waterplaneMoments(testCase.MeshCube.ConnectivityList,testCase.MeshCube.Points,FN);

            testCase.verifyLessThan(abs(Sy),testCase.epsilon);
        end

        function testWaterplaneMoments4(testCase)
            FN = -faceNormal(testCase.MeshCube);
            [~,~,~,Sxy,~,~] = waterplaneMoments(testCase.MeshCube.ConnectivityList,testCase.MeshCube.Points,FN);

            testCase.verifyLessThan(abs(Sxy),testCase.epsilon);
        end

        function testWaterplaneMoments5(testCase)
            FN = -faceNormal(testCase.MeshCube);
            [~,~,~,~,Sxx,~] = waterplaneMoments(testCase.MeshCube.ConnectivityList,testCase.MeshCube.Points,FN);

            testCase.verifyLessThan(abs(Sxx),testCase.epsilon);
        end

        function testWaterplaneMoments6(testCase)
            FN = -faceNormal(testCase.MeshCube);
            [~,~,~,~,~,Syy] = waterplaneMoments(testCase.MeshCube.ConnectivityList,testCase.MeshCube.Points,FN);

            testCase.verifyLessThan(abs(Syy),testCase.epsilon);
        end

        function testWaterplaneMoments7(testCase)
            FN = faceNormal(testCase.ModelSquare);
            [S,~,~,~,~,~] = waterplaneMoments(testCase.ModelSquare.ConnectivityList,testCase.ModelSquare.Points,FN);
            value = 1;

            testCase.verifyLessThan(abs(S - value),testCase.epsilon);
        end

        function testWaterplaneMoments8(testCase)
            FN = faceNormal(testCase.ModelSquare);
            [~,Sx,~,~,~,~] = waterplaneMoments(testCase.ModelSquare.ConnectivityList,testCase.ModelSquare.Points,FN);
            value = 0.5;

            testCase.verifyLessThan(abs(Sx - value),testCase.epsilon);
        end

        function testWaterplaneMoments9(testCase)
            FN = faceNormal(testCase.ModelSquare);
            [~,~,Sy,~,~,~] = waterplaneMoments(testCase.ModelSquare.ConnectivityList,testCase.ModelSquare.Points,FN);
            value = 0.5;

            testCase.verifyLessThan(abs(Sy - value),testCase.epsilon);
        end

        function testWaterplaneMoments10(testCase)
            FN = faceNormal(testCase.ModelSquare);
            [~,~,~,Sxy,~,~] = waterplaneMoments(testCase.ModelSquare.ConnectivityList,testCase.ModelSquare.Points,FN);
            value = 0.25;

            testCase.verifyLessThan(abs(Sxy - value),testCase.epsilon);
        end

        function testWaterplaneMoments11(testCase)
            FN = faceNormal(testCase.ModelSquare);
            [~,~,~,~,Sxx,~] = waterplaneMoments(testCase.ModelSquare.ConnectivityList,testCase.ModelSquare.Points,FN);
            value = 1/3;

            testCase.verifyLessThan(abs(Sxx - value),testCase.epsilon);
        end

        function testWaterplaneMoments12(testCase)
            FN = faceNormal(testCase.ModelSquare);
            [~,~,~,~,~,Syy] = waterplaneMoments(testCase.ModelSquare.ConnectivityList,testCase.ModelSquare.Points,FN);
            value = 1/3;

            testCase.verifyLessThan(abs(Syy - value),testCase.epsilon);
        end

        function testWaterplaneMoments13(testCase)
            FN = faceNormal(testCase.ModelSquare2);
            [S,~,~,~,~,~] = waterplaneMoments(testCase.ModelSquare2.ConnectivityList,testCase.ModelSquare2.Points,FN);
            value = 9;

            testCase.verifyLessThan(abs(S - value),testCase.epsilon);
        end

        function testWaterplaneMoments14(testCase)
            FN = faceNormal(testCase.ModelSquare2);
            [~,Sx,~,~,~,~] = waterplaneMoments(testCase.ModelSquare2.ConnectivityList,testCase.ModelSquare2.Points,FN);
            value = 22.5;

            testCase.verifyLessThan(abs(Sx - value),testCase.epsilon);
        end

        function testWaterplaneMoments15(testCase)
            FN = faceNormal(testCase.ModelSquare2);
            [~,~,Sy,~,~,~] = waterplaneMoments(testCase.ModelSquare2.ConnectivityList,testCase.ModelSquare2.Points,FN);
            value = 27;

            testCase.verifyLessThan(abs(Sy - value),testCase.epsilon);
        end

        function testWaterplaneMoments16(testCase)
            FN = faceNormal(testCase.ModelSquare2);
            [~,~,~,Sxy,~,~] = waterplaneMoments(testCase.ModelSquare2.ConnectivityList,testCase.ModelSquare2.Points,FN);
            value = 69.75;

            testCase.verifyLessThan(abs(Sxy - value),testCase.epsilon);
        end

        function testWaterplaneMoments17(testCase)
            FN = faceNormal(testCase.ModelSquare2);
            [~,~,~,~,Sxx,~] = waterplaneMoments(testCase.ModelSquare2.ConnectivityList,testCase.ModelSquare2.Points,FN);
            value = 63;

            testCase.verifyLessThan(abs(Sxx - value),testCase.epsilon);
        end

        function testWaterplaneMoments18(testCase)
            FN = faceNormal(testCase.ModelSquare2);
            [~,~,~,~,~,Syy] = waterplaneMoments(testCase.ModelSquare2.ConnectivityList,testCase.ModelSquare2.Points,FN);
            value = 88.5;

            testCase.verifyLessThan(abs(Syy - value),testCase.epsilon);
        end
    end
end

