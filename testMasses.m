classdef testMasses < matlab.unittest.TestCase

    properties
        epsilon
    end    
 
    methods(TestMethodSetup)
        function defineEpsilon(testCase)
            testCase.epsilon = 1E-6;
        end
    end
 
    methods(TestMethodTeardown)
    end
    
    methods (Test)
        function testMassMoments1(testCase)
            [m,~] = massMoments( ...
                [1,randn(1,3); ...
                2,randn(1,3)]);
            value = 3;

            testCase.verifyEqual(m,value);
        end

        function testMassMoments2(testCase)
            [m,~] = massMoments( ...
                [-2,randn(1,3); ...
                2,randn(1,3)]);
            value = 0;

            testCase.verifyEqual(m,value);
        end

        function testMassMoments3(testCase)
            [m,~] = massMoments( ...
                [1,randn(1,3); ...
                2,randn(1,3); ...
                20,randn(1,3); ...
                2,randn(1,3); ...
                -1,randn(1,3); ...
                2.2,randn(1,3)]);
            value = 26.2;

            testCase.verifyLessThan(abs(m - value), testCase.epsilon);
        end

        function testMassMoments4(testCase)
            [~,Cg] = massMoments( ...
                [1,0,0,0; ...
                3,0,0,1;]);
            value = [0,0,0.75];

            testCase.verifyLessThan(abs(Cg - value), testCase.epsilon);
        end

        function testMassMoments5(testCase)
            [~,Cg] = massMoments( ...
                [9,-1,3,-3; ...
                9,3,1,1;]);
            value = [1,2,-1];

            testCase.verifyLessThan(abs(Cg - value), testCase.epsilon);
        end

        function testMassMoments6(testCase)
            [~,Cg] = massMoments( ...
                [2,10,10,10; ...
                1,5,5,5; ...
                1,15,15,15;]);
            value = [10,10,10];

            testCase.verifyLessThan(abs(Cg - value), testCase.epsilon);
        end
    end
end