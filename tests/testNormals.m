classdef testNormals < matlab.unittest.TestCase

    properties
        epsilon
    end    
 
    methods(TestMethodSetup)
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
        function testNormal6DOF1(testCase)
            [FN6] = normal6DOF([1,1,1],[0,0,0]);
            value = [0,0,0,0,0,0];

            testCase.verifyEqual(FN6,value);
        end

        function testNormal6DOF2(testCase)
            [FN6] = normal6DOF([0,0,0],[1,1,1]);
            value = [1,1,1,0,0,0];

            testCase.verifyEqual(FN6,value);
        end

        function testNormal6DOF3(testCase)
            [FN6] = normal6DOF([1,1,1],[2,3,4]);
            value = [2,3,4,1,-2,1];

            testCase.verifyEqual(FN6,value);
        end

        function testNormal6DOF4(testCase)
            [FN6] = normal6DOF([2,3,4],[1,1,1]);
            value = [1,1,1,-1,2,-1];

            testCase.verifyEqual(FN6,value);
        end
    end
end