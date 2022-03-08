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
            w = 0;
            p = 1026;
            FN = [1,1,1] / norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,1,0;
                   0,0,1];
            [~,B] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = 0;

            testCase.verifyLessThan(abs(B - value),testCase.epsilon)
        end

        function testAddedMassAndDamping2(testCase)
            w = 0;
            p = 1026;
            FN = [1,1,1] / norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,1,0;
                   0,0,1];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = FN(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping3(testCase)
            w = 0;
            p = 1026;
            FN = [1,1,1] / norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [0,1,0;
                   1,0,0;
                   0,0,1];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = FN(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping4(testCase)
            w = 0;
            p = 1026;
            FN = [1,1,1] / norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = FN(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping5(testCase)
            w = 6;
            p = 1026;
            FN = [1,1,1] / norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = FN(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping6(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [1;1;1;1;1;1];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = FN(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping7(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [5;3;7;2;8;9]; %[1+7i;7;3i;4;5;-6];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = FN(i)*phi(j)*p*0.8660254037844389;
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping8(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [5 + 4i;3 + 2i;7 + 0.1i;2 + 6i;8 + 3.3i;9]; %[1+7i;7;3i;4;5;-6];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = real(FN(i)*phi(j)*p*0.8660254037844389);
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping9(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [1+7i;7;3i;4;5;-6];
            Tri = [0,1,0;
                   0,0,0;
                   1,0,0];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = real(FN(i)*phi(j)*p*0.5);
                end
            end
            
            testCase.verifyLessThan(abs(A(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping10(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [1+7i;7;3i;4;5;-6];
            Tri = [0,1,0;
                   0,0,0;
                   1,0,0];
            [A,~] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = real(FN(i)*phi(j)*p*0.5);
                end

                value(4,j) = real(FN(3)/6*phi(j)*p);
                value(5,j) = real(-FN(3)/6*phi(j)*p);
                value(6,j) = real((FN(2)-FN(1))/6*phi(j)*p);
            end
            
            testCase.verifyLessThan(abs(A(:,:) - value(:,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping11(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [1+7i;7;3i;4;5;-6];
            Tri = [0,1,0;
                   0,0,0;
                   1,0,0];
            [~,B] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = imag(FN(i)*phi(j)*p*0.5);
                end

                value(4,j) = imag(FN(3)/6*phi(j)*p);
                value(5,j) = imag(-FN(3)/6*phi(j)*p);
                value(6,j) = imag((FN(2)-FN(1))/6*phi(j)*p);
            end

            value = -value*w; 
            
            testCase.verifyLessThan(abs(B - value),testCase.epsilon)
        end

        function testAddedMassAndDamping12(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [1+7i;7;3i;4;5;-6];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            [~,B] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = imag(FN(i)*phi(j)*p*0.8660254037844389);
                end
            end

            value = -value*w; 
            
            testCase.verifyLessThan(abs(B(1:3,:) - value(1:3,:)),testCase.epsilon)
        end

        function testAddedMassAndDamping13(testCase)
            w = 6;
            p = 1026;
            FN = [2,3,5] / norm([2,3,5]);
            phi = [1-7i;7;-3i;-4 + 1i;5;-6];
            Tri = [0,1,0;
                   0,0,1;
                   1,0,0];
            [~,B] = addedMassAndDampingMatrices(Tri,phi,FN,p,w);
            value = zeros(6);
            for j = 1:6
                for i = 1:3
                    value(i,j) = imag(FN(i)*phi(j)*p*0.8660254037844389);
                end
            end

            value = -value*w; 
            
            testCase.verifyLessThan(abs(B(1:3,:) - value(1:3,:)),testCase.epsilon)
        end
    end
end

