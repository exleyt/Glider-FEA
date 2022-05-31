classdef testAddedMassAndDamping < matlab.unittest.TestCase

    properties
        epsilon
    end

    methods(TestMethodSetup)
        function defineEpsilon(testCase)
            testCase.epsilon = 1E-4;
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

        function testAddedMassAndDamping1(testCase)
            p = 1026;
            FN = [1,1,1] / norm([1,1,1]);
            FN6 = normal6DOF([0,0,0],FN);
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,1,0;
                   0,0,1];
            A = addedMassAndDampingMatrices(Tri,phi,FN6,p);
            value = zeros(6);
            for j = 1:6
                for i = 1:6
                    value(i,j) = FN6(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A - value),testCase.epsilon)
        end

        function testAddedMassAndDamping2(testCase)
            p = 1026;
            FN = [1,1,1] / norm([1,1,1]);
            FN6 = normal6DOF([1,0.5,-1],FN);
            phi = [1;1;1;1;1;1];
            Tri = [0,1,0;
                   1,0,0;
                   0,0,1];
            A = addedMassAndDampingMatrices(Tri,phi,FN6,p);
            value = zeros(6);
            for j = 1:6
                for i = 1:6
                    value(i,j) = FN6(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A - value),testCase.epsilon)
        end

        function testAddedMassAndDamping3(testCase)
            p = 1020;
            FN = [1,1,1] / norm([1,1,1]);
            FN6 = normal6DOF([1,1,-1.1],FN);
            phi = [1;1;1;1;1;1];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            A = addedMassAndDampingMatrices(Tri,phi,FN6,p);
            value = zeros(6);
            for j = 1:6
                for i = 1:6
                    value(i,j) = FN6(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A - value),testCase.epsilon)
        end

        function testAddedMassAndDamping4(testCase)
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            FN6 = normal6DOF([0.2,-0.3,-0.1],FN);
            phi = [1;1;1;1;1;1];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            A = addedMassAndDampingMatrices(Tri,phi,FN6,p);
            value = zeros(6);
            for j = 1:6
                for i = 1:6
                    value(i,j) = FN6(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A - value),testCase.epsilon)
        end

        function testAddedMassAndDamping5(testCase)
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            FN6 = normal6DOF([0.2,-0.3,-0.1],FN);
            phi = [5;3;7;2;8;9];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            A = addedMassAndDampingMatrices(Tri,phi,FN6,p);
            value = zeros(6);
            for j = 1:6
                for i = 1:6
                    value(i,j) = FN6(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A - value),testCase.epsilon)
        end

        function testAddedMassAndDamping6(testCase)
            p = 1021;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [5 + 4i;3 + 2i;7 + 0.1i;2 + 6i;8 + 3.3i;9];
            FN6 = normal6DOF([0.2,-0.3,-0.1],FN);
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            A = addedMassAndDampingMatrices(Tri,phi,FN6,p);
            value = zeros(6);
            for j = 1:6
                for i = 1:6
                    value(i,j) = FN6(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A - value),testCase.epsilon)
        end

        function testAddedMassAndDamping7(testCase)
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            FN6 = normal6DOF([0.2,-0.3,-0.1],FN);
            phi = [1+7i;7;3i;4;5;-6];
            Tri = [0,1,0;
                   0,0,0;
                   1,0,0];
            A = addedMassAndDampingMatrices(Tri,phi,FN6,p);
            value = zeros(6);
            for j = 1:6
                for i = 1:6
                    value(i,j) = FN6(i)*phi(j)*p*0.5;
                end
            end
            
            testCase.verifyLessThan(abs(A - value),testCase.epsilon)
        end
end

