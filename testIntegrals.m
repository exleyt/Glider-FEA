classdef testIntegrals < matlab.unittest.TestCase

    properties
        epsilon
        helper
    end    
 
    methods(TestMethodSetup)
        function defineEpsilon(testCase)
            testCase.epsilon = 1E-4;
        end

        function defineHelper(testCase)
            testCase.helper = testIntegralsHelper;
        end
    end
 
    methods(TestMethodTeardown)
    end
    
    methods (Test)
        function testSurfIntPhiI1(testCase)
            su = 1;
            sv = 1;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI2(testCase)
            su = 0;
            sv = 10;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI3(testCase)
            su = 0;
            sv = -10;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI4(testCase)
            su = 6;
            sv = 0;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI5(testCase)
            su = -6;
            sv = 0;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI6(testCase)
            su = 0;
            sv = 0;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI7(testCase)
            su = -1.5;
            sv = -90;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI8(testCase)
            su = -1.5;
            sv = 7;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI9(testCase)
            su = 2;
            sv = -13;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end

        function testSurfIntPhiI10(testCase)
            su = 15;
            sv = 3;

            [a1,au,av] = surfIntPhiI(su,sv);
            [b1,bu,bv] = testCase.helper.evaluateIntegral(su,sv);

            testCase.verifyLessThan(abs(a1-b1),testCase.epsilon);
            testCase.verifyLessThan(abs(au-bu),testCase.epsilon);
            testCase.verifyLessThan(abs(av-bv),testCase.epsilon);
        end
    end
end