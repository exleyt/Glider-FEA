classdef testASQ < matlab.unittest.TestCase
    
    properties
        epsilon
        epsilonASQ
        helper
    end    
 
    methods(TestMethodSetup)
        function defineEpsilon(testCase)
            testCase.epsilon = 1E-6;
            testCase.epsilonASQ = 1E-2;
        end

        function defineHelper(testCase)
            testCase.helper = testIntegralsHelper;
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
        function testPASQ1(testCase)
            f = @(x) x;
            a = 0;
            b = 1;
            int = integral(f,a,b);
            value = testCase.helper.pasq(f,a,b);
            
            testCase.verifyLessThan(abs(int - value), testCase.epsilon)
        end

        function testPASQ2(testCase)
            f = @(x) x.^2;
            a = -1;
            b = 1;
            int = integral(f,a,b);
            value = testCase.helper.pasq(f,a,b);
            
            testCase.verifyLessThan(abs(int - value), testCase.epsilon)
        end

        function testPASQ3(testCase)
            f = @(x) x.^3;
            a = -2;
            b = 1;
            int = integral(f,a,b);
            value = testCase.helper.pasq(f,a,b);
            
            testCase.verifyLessThan(abs(int - value), testCase.epsilon)
        end

        function testPASQ4(testCase)
            f = @(x) x + x.^2 + x.^3;
            a = -2.5;
            b = 0;
            int = integral(f,a,b);
            value = testCase.helper.pasq(f,a,b);
            
            testCase.verifyLessThan(abs(int - value), testCase.epsilon)
        end

        function testPASQ5(testCase)
            f = @(x) x.*x.^2 - x.^3;
            a = 1.1;
            b = 1.3;
            int = integral(f,a,b);
            value = testCase.helper.pasq(f,a,b);
            
            testCase.verifyLessThan(abs(int - value), testCase.epsilon)
        end

        function testASQ1(testCase)
            R = 0;
            sx3 = 0;
            K = 0;
            f = @(x) besselj(0,x*R)*exp(x*sx3)/(x-K);
            a = 1.1;
            b = 1.3;
            pasq = testCase.helper.pasq(f,a,b);
            value = asq(R,sx3,K,a,b);           
            testCase.verifyLessThan(abs(pasq - value), testCase.epsilonASQ)
        end

        function testASQ2(testCase)
            R = 0;
            sx3 = 0;
            K = 1;
            f = @(x) besselj(0,x*R)*exp(x*sx3)/(x-K);
            a = 1.1;
            b = 1.3;
            pasq = testCase.helper.pasq(f,a,b);
            value = asq(R,sx3,K,a,b);           
            testCase.verifyLessThan(abs(pasq - value), testCase.epsilonASQ)
        end

        function testASQ3(testCase)
            R = 0;
            sx3 = -2;
            K = 1;
            f = @(x) besselj(0,x*R)*exp(x*sx3)/(x-K);
            a = 1.1;
            b = 1.3;
            pasq = testCase.helper.pasq(f,a,b);
            value = asq(R,sx3,K,a,b);           
            testCase.verifyLessThan(abs(pasq - value), testCase.epsilonASQ)
        end

        function testASQ4(testCase)
            R = 1.2;
            sx3 = -2;
            K = 1;
            f = @(x) besselj(0,x*R)*exp(x*sx3)/(x-K);
            a = 1.1;
            b = 1.3;
            pasq = testCase.helper.pasq(f,a,b);
            value = asq(R,sx3,K,a,b);           
            testCase.verifyLessThan(abs(pasq - value), testCase.epsilonASQ)
        end

        function testASQ5(testCase)
            R = 0.2;
            sx3 = -1.1;
            K = 4;
            f = @(x) besselj(0,x*R)*exp(x*sx3)/(x-K);
            a = 1.2;
            b = 3;
            pasq = testCase.helper.pasq(f,a,b);
            value = asq(R,sx3,K,a,b);           
            testCase.verifyLessThan(abs(pasq - value), testCase.epsilonASQ)
        end
    end
end