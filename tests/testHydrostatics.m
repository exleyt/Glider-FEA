classdef testHydrostatics < matlab.unittest.TestCase

    properties
        ModelSqr1
        ModelSqr1Off1
        ModelSqr1Off231
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

        function createModelSqr1Off1(testCase)
            testCase.ModelSqr1Off1 = createpde();
            W = 1;
            D = 1;
            H = 1;
            testCase.ModelSqr1Off1.Geometry = multicuboid(W,D,H);
            translate(testCase.ModelSqr1Off1.Geometry,[1,1,1]);
            generateMesh(testCase.ModelSqr1Off1,'GeometricOrder','linear');
        end

        function createModelSqr1Off231(testCase)
            testCase.ModelSqr1Off231 = createpde();
            W = 1;
            D = 1;
            H = 1;
            testCase.ModelSqr1Off231.Geometry = multicuboid(W,D,H);
            translate(testCase.ModelSqr1Off231.Geometry,[2,3,1]);
            generateMesh(testCase.ModelSqr1Off231,'GeometricOrder','linear');
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
        function testHydrostaticRestoringMatrix(testCase)
            pm = [10,0,0,0];
            p = 1026;
            g = 9.81;
            model = testCase.ModelSqr1;
            face = 2;
            to = triangulation(model.Mesh.Elements.',model.Mesh.Nodes.');
            [T,P] = freeBoundary(to);
            fto = triangulation(T,P);
            CP = incenter(fto);
            FaceIDs = nearestFace(model.Geometry,CP(:,:));
            [~,N] = size(FaceIDs);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                end
            end
            CL = zeros(c,3);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                    CL(c,:) = T(i,:);
                end
            end
            mesh = triangulation(CL,P);
            CB = [0,0,0];
            V = 1;
            C = hydrostaticRestoringMatrix(pm,g,p,V,CB,mesh);
            Ch = p*g*1;
            c44 = p*g*1/12;
            c55 = p*g*1/12;
            value = [
                0,  0,  0,          0,                  0,                  0;
                0,  0,  0,          0,                  0,                  0;
                0,  0,  -Ch,         0,                  0,                  0;
                0,  0,  0,          -c44,                0,                  0;
                0,  0,  0,          0,                  -c55,                0;
                0,  0,  0,          0,                  0,                  0;
            ];
            
            testCase.verifyLessThan(abs(C - value), testCase.epsilon)
        end

        function testHydrostaticRestoringMatrix2(testCase)
            pm = [10,3,2,1];
            p = 1026;
            g = 9.81;
            model = testCase.ModelSqr1;
            face = 2;
            to = triangulation(model.Mesh.Elements.',model.Mesh.Nodes.');
            [T,P] = freeBoundary(to);
            fto = triangulation(T,P);
            CP = incenter(fto);
            FaceIDs = nearestFace(model.Geometry,CP(:,:));
            [~,N] = size(FaceIDs);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                end
            end
            CL = zeros(c,3);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                    CL(c,:) = T(i,:);
                end
            end
            mesh = triangulation(CL,P);
            CB = [0,0,0];
            V = 1;
            C = hydrostaticRestoringMatrix(pm,g,p,V,CB,mesh);
            Ch = p*g*1;
            W = 10*g;
            c44 = -W*1 - p*g*1/12;
            c55 = -W*1 - p*g*1/12;
            value = [
                0,  0,  0,          0,                  0,                  0;
                0,  0,  0,          0,                  0,                  0;
                0,  0,  -Ch,         0,                  0,                 0;
                0,  0,  0,          c44,                0,                  W;
                0,  0,  0,          0,                  c55,                W*2;
                0,  0,  0,          0,                  0,                  0;
            ];
            
            testCase.verifyLessThan(abs(C - value), testCase.epsilon)
        end

        function testHydrostaticRestoringMatrix3(testCase)
            pm = [11,2,3,2];
            p = 1026;
            g = 9.81;
            model = testCase.ModelRect243;
            face = 2;
            to = triangulation(model.Mesh.Elements.',model.Mesh.Nodes.');
            [T,P] = freeBoundary(to);
            fto = triangulation(T,P);
            CP = incenter(fto);
            FaceIDs = nearestFace(model.Geometry,CP(:,:));
            [~,N] = size(FaceIDs);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                end
            end
            CL = zeros(c,3);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                    CL(c,:) = T(i,:);
                end
            end
            mesh = triangulation(CL,P);
            CB = [1,2,1.5];
            V = 24;
            C = hydrostaticRestoringMatrix(pm,g,p,V,CB,mesh);
            Ch = p*g*8;
            Fb =  p*g*24;
            W = 11*g;
            c44 = -W*2 + Fb*CB(3) - p*g*32/3;
            c55 = -W*2 + Fb*CB(3) - p*g*8/3;
            value = [
                0,  0,  0,          0,                  0,                  0;
                0,  0,  0,          0,                  0,                  0;
                0,  0,  -Ch,         0,                  0,                 0;
                0,  0,  0,          c44,                0,                  -Fb*CB(1) + W*2;
                0,  0,  0,          0,                  c55,                -Fb*CB(2) + W*3;
                0,  0,  0,          0,                  0,                  0;
            ];
            
            testCase.verifyLessThan(abs(C - value), testCase.epsilon)
        end

        function testHydrostaticRestoringMatrix4(testCase)
            pm = [10,-1,2,1;
                  10,1,-2,-1];
            p = 1026;
            g = 9.81;
            model = testCase.ModelRect243;
            face = 2;
            to = triangulation(model.Mesh.Elements.',model.Mesh.Nodes.');
            [T,P] = freeBoundary(to);
            fto = triangulation(T,P);
            CP = incenter(fto);
            FaceIDs = nearestFace(model.Geometry,CP(:,:));
            [~,N] = size(FaceIDs);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                end
            end
            CL = zeros(c,3);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                    CL(c,:) = T(i,:);
                end
            end
            mesh = triangulation(CL,P);
            CB = [1,2,1.5];
            V = 24;
            C = hydrostaticRestoringMatrix(pm,g,p,V,CB,mesh);
            Ch = p*g*8;
            Fb =  p*g*24;
            c44 = 0 + Fb*CB(3) - p*g*32/3;
            c55 = 0 + Fb*CB(3) - p*g*8/3;
            value = [
                0,  0,  0,          0,                  0,                  0;
                0,  0,  0,          0,                  0,                  0;
                0,  0,  -Ch,         0,                  0,                  0;
                0,  0,  0,          c44,                0,                  -Fb*CB(1);
                0,  0,  0,          0,                  c55,                -Fb*CB(2);
                0,  0,  0,          0,                  0,                  0;
            ];
            
            testCase.verifyLessThan(abs(C - value), testCase.epsilon)
        end

        function testHydrostaticRestoringMatrix5(testCase)
            pm = [2,3,4,5];
            p = 1026;
            g = 9.81;
            model = testCase.ModelSqr1Off1;
            face = 2;
            to = triangulation(model.Mesh.Elements.',model.Mesh.Nodes.');
            [T,P] = freeBoundary(to);
            fto = triangulation(T,P);
            CP = incenter(fto);
            FaceIDs = nearestFace(model.Geometry,CP(:,:));
            [~,N] = size(FaceIDs);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                end
            end
            CL = zeros(c,3);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                    CL(c,:) = T(i,:);
                end
            end
            mesh = triangulation(CL,P);
            CB = [1,1,1];
            V = 1;
            C = hydrostaticRestoringMatrix(pm,g,p,V,CB,mesh);
            Ch = p*g*1;
            Fb =  p*g*1;
            W = 2*g;
            c44 = -W*5 + Fb*1 - p*g*(1 + 1/12);
            c55 = -W*5 + Fb*1 - p*g*(1 + 1/12);
            value = [
                0,  0,  0,          0,                  0,                  0;
                0,  0,  0,          0,                  0,                  0;
                0,  0,  -Ch,         -Ch*1,               Ch*1,              0;
                0,  0,  -Ch*1,       c44,                Ch*1*1,            W*5 - Fb*1;
                0,  0,  Ch*1,      Ch*1*1,            c55,                W*4 - Fb*1;
                0,  0,  0,          0,                  0,                  0;
            ];
            
            testCase.verifyLessThan(abs(C - value), testCase.epsilon)
        end

        function testHydrostaticRestoringMatrix6(testCase)
            pm = [2,3,4,5];
            p = 1026;
            g = 9.81;
            model = testCase.ModelSqr1Off231;
            face = 2;
                        to = triangulation(model.Mesh.Elements.',model.Mesh.Nodes.');
            [T,P] = freeBoundary(to);
            fto = triangulation(T,P);
            CP = incenter(fto);
            FaceIDs = nearestFace(model.Geometry,CP(:,:));
            [~,N] = size(FaceIDs);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                end
            end
            CL = zeros(c,3);
            c = 0;
            for i = 1:N
                if FaceIDs(1,i) == face
                    c = c + 1;
                    CL(c,:) = T(i,:);
                end
            end
            mesh = triangulation(CL,P);
            CB = [2,3,2.5];
            V = 6;
            C = hydrostaticRestoringMatrix(pm,g,p,V,CB,mesh);
            Ch = p*g*1;
            Fb =  p*g*6;
            W = 2*g;
            c44 = -W*5 + Fb*2.5 - p*g*(9 + 1/12);
            c55 = -W*5 + Fb*2.5 - p*g*(4 + 1/12);
            value = [
                0,  0,  0,          0,                  0,                  0;
                0,  0,  0,          0,                  0,                  0;
                0,  0,  -Ch,         -Ch*3,               Ch*2,              0;
                0,  0,  -Ch*3,       c44,                Ch*2*3,            W*5 - Fb*2;
                0,  0,  Ch*2,      Ch*2*3,            c55,                W*4 - Fb*3;
                0,  0,  0,          0,                  0,                  0;
            ];
            
            testCase.verifyLessThan(abs(C - value), testCase.epsilon)
        end
    end
end

