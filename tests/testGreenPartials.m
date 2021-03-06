classdef testGreenPartials < matlab.unittest.TestCase

    properties
        epsilon
    end    
 
    methods(TestMethodSetup)
        function defineEpsilon(testCase)
            testCase.epsilon = 1E-2;
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
        function testGreenFunctionAndPartialXINormal1(testCase)
            x = [-1,8,-1];
            xi = [-6,-5.4,-2];
            N = [3,3,3];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*testCase.epsilon^2,K);
            [f,~] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = G(0);
            
            testCase.verifyLessThan(abs(f - value), testCase.epsilon)
        end

        function testGreenFunctionAndPartialXINormal2(testCase)
            x = [.1,.8,-.1];
            xi = [.6,-5.4,-.2];
            N = [3,3,3];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*testCase.epsilon^2,K);
            [f,~] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = G(0);
            
            testCase.verifyLessThan(abs(f - value), testCase.epsilon)
        end

        function testGreenFunctionAndPartialXINormal3(testCase)
            x = [.1,.8,-.1];
            xi = [.6,-5.4,-.2];
            N = [3,3,3];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*y,K);
            [~,df] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = (G(testCase.epsilon^2) - G(0))/testCase.epsilon^2;
            
            testCase.verifyLessThan(abs(df - value), testCase.epsilon*(1+1i))
        end

        function testGreenFunctionAndPartialXINormal4(testCase)
            x = [1.1,2.8,-4.1];
            xi = [2.6,-0.4,-3.2];
            N = [3,3,3];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*y,K);
            [~,df] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = (G(testCase.epsilon^2) - G(0))/testCase.epsilon^2;
            
            testCase.verifyLessThan(abs(df - value), testCase.epsilon*(1+1i))
        end

        function testGreenFunctionAndPartialXINormal5(testCase)
            x = [.7,-6.2,-1];
            xi = [.6,-5.4,-.2];
            N = [3,3,3];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*y,K);
            [~,df] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = (G(testCase.epsilon^2) - G(0))/testCase.epsilon^2;
            
            testCase.verifyLessThan(abs(df - value), testCase.epsilon*(1+1i))
        end

        function testGreenFunctionAndPartialXINormal6(testCase)
            x = [-1,8,-1];
            xi = [-6,-5.4,-2];
            N = [2,3,4];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*testCase.epsilon^2,K);
            [f,~] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = G(0);
            
            testCase.verifyLessThan(abs(f - value), testCase.epsilon)
        end

        function testGreenFunctionAndPartialXINormal7(testCase)
            x = [.1,.8,-.1];
            xi = [.6,-5.4,-.2];
            N = [9,4,2];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*testCase.epsilon^2,K);
            [f,~] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = G(0);
            
            testCase.verifyLessThan(abs(f - value), testCase.epsilon)
        end

        function testGreenFunctionAndPartialXINormal8(testCase)
            x = [.1,.8,-.1];
            xi = [.6,-5.4,-.2];
            N = [-1,-2,-3];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*y,K);
            [~,df] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = (G(testCase.epsilon^2) - G(0))/testCase.epsilon^2;
            
            testCase.verifyLessThan(abs(df - value), testCase.epsilon*(1+1i)*10)
        end

        function testGreenFunctionAndPartialXINormal9(testCase)
            x = [1.1,2.8,-4.1];
            xi = [2.6,-0.4,-3.2];
            N = [-3,1,-5];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*y,K);
            [~,df] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = (G(testCase.epsilon^2) - G(0))/testCase.epsilon^2;
            
            testCase.verifyLessThan(abs(df - value), testCase.epsilon*(1+1i))
        end

        function testGreenFunctionAndPartialXINormal10(testCase)
            x = [.7,-6.2,-1];
            xi = [.6,-5.4,-.2];
            N = [3,-3.5,-3.2];
            K = 1.2;
            G = @(y) greenFunction(x,xi + N*y,K);
            [~,df] = greenFunctionAndPartialXINormal(x,xi,N,K);
            value = (G(testCase.epsilon^2) - G(0))/testCase.epsilon^2;
            
            testCase.verifyLessThan(abs(df - value), testCase.epsilon*(1+1i)*10)
        end
    end
end

