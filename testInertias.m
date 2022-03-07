classdef testInertias < matlab.unittest.TestCase

    properties
        epsilon
    end

    methods(TestMethodSetup)
        function defineEpsilon(testCase)
            testCase.epsilon = 1E-4;
        end
    end
    
    methods (Test)
        function testMomentOfInertia1(testCase)
            pm = [2,1,0,0];
            I = momentOfInertia(pm);
            value = [0,0,0;
                    0,2,0;
                    0,0,2];

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia2(testCase)
            pm = [4,1,0,0;
                  2,3,0,0];
            I = momentOfInertia(pm);
            value = [0,0,0;
                    0,2*3^2 + 4*1^2,0;
                    0,0,2*3^2 + 4*1^2];

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia3(testCase)
            pm = [2,0,1,0];
            I = momentOfInertia(pm);
            value = [2,0,0;
                    0,0,0;
                    0,0,2];

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia4(testCase)
            pm = [2,0,3,0;
                  2,0,1,0];
            I = momentOfInertia(pm);
            value = [2*3^2 + 2*1^2,0,0;
                    0,0,0;
                    0,0,2*3^2 + 2*1^2];

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia5(testCase)
            pm = [7,0,0,1];
            I = momentOfInertia(pm);
            value = [7,0,0;
                    0,7,0;
                    0,0,0];

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia6(testCase)
            pm = [2,0,0,3;
                  5,0,0,1];
            I = momentOfInertia(pm);
            value = [2*3^2 + 5*1^2,0,0;
                    0,2*3^2 + 5*1^2,0;
                    0,0,0];

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia7(testCase)
            pm = [3,1,2,0];
            I = momentOfInertia(pm);
            value = [3*2^2,-3*1*2,0;
                    -3*1*2,3*1^2,0;
                    0,0,3*(1^2 + 2^2)];

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia8(testCase)
            x = 1.7;
            y = -5.5;
            z = -13.2;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = momentOfInertia(pm);
            value = [y^2 + z^2,-x*y,-x*z;
                    -x*y,x^2 + z^2,-y*z;
                    -x*z,-y*z,x^2 + y^2]*m;

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testMomentOfInertia9(testCase)
            x = -1.17;
            y = 9.0;
            z = 0.1;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = momentOfInertia(pm);
            value = [y^2 + z^2,-x*y,-x*z;
                    -x*y,x^2 + z^2,-y*z;
                    -x*z,-y*z,x^2 + y^2]*m;

            testCase.verifyLessThan(abs(I - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix1(testCase)
            x = -1.17;
            y = 9.0;
            z = 0.1;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1 + z - y;
                    1 - z + x;
                    1 + y - x]*m;

            testCase.verifyLessThan(abs(sum(I(1:3,:),2) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix2(testCase)
            x = 7;
            y = 0;
            z = 0;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1 + z - y;
                    1 - z + x;
                    1 + y - x]*m;

            testCase.verifyLessThan(abs(sum(I(1:3,:),2) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix3(testCase)
            x = 0;
            y = 3;
            z = 0;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1 + z - y;
                    1 - z + x;
                    1 + y - x]*m;

            testCase.verifyLessThan(abs(sum(I(1:3,:),2) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix4(testCase)
            x = 0;
            y = 0;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1 + z - y;
                    1 - z + x;
                    1 + y - x]*m;

            testCase.verifyLessThan(abs(sum(I(1:3,:),2) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix5(testCase)
            x = 2;
            y = 2;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1 + z - y;
                    1 - z + x;
                    1 + y - x]*m;

            testCase.verifyLessThan(abs(sum(I(1:3,:),2) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix6(testCase)
            x = 3;
            y = 2;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1;
                m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1 + z - y;
                    1 - z + x;
                    1 + y - x]*m*2;

            testCase.verifyLessThan(abs(sum(I(1:3,:),2) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix7(testCase)
            x = 3;
            y = 2;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,-R1;
                m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1 + 0;
                    1 - 0;
                    1 + 0]*m*2;

            testCase.verifyLessThan(abs(sum(I(1:3,:),2) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix9(testCase)
            x = 30;
            y = 20;
            z = -40;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1;
                m,R1];
            I = bodyInertiaMatrix(pm);
            value = [1,0,0,0,z,-y;
                    0,1,0,-z,0,x;
                    0,0,1,y,-x,0]*m*2;

            testCase.verifyLessThan(abs(I(1:3,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix10(testCase)
            x = -1.17;
            y = 9.0;
            z = 0.1;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,-z,y;
                    z,0,-x;
                    -y,x,0]*m,momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix11(testCase)
            x = 7;
            y = 0;
            z = 0;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,-z,y;
                    z,0,-x;
                    -y,x,0]*m,momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix12(testCase)
            x = 0;
            y = 3;
            z = 0;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,-z,y;
                    z,0,-x;
                    -y,x,0]*m,momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix13(testCase)
            x = 0;
            y = 0;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,-z,y;
                    z,0,-x;
                    -y,x,0]*m,momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix14(testCase)
            x = 2;
            y = 2;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,-z,y;
                    z,0,-x;
                    -y,x,0]*m,momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix15(testCase)
            x = 3;
            y = 2;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1;
                m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,-z,y;
                    z,0,-x;
                    -y,x,0]*m*2,momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix16(testCase)
            x = 3;
            y = 2;
            z = -4;
            R1 = [x,y,z];
            m = 3;
            pm = [m,-R1;
                m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,0,0;
                    0,0,0;
                    0,0,0],momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end

        function testBodyInertiaMatrix17(testCase)
            x = 30;
            y = 20;
            z = -40;
            R1 = [x,y,z];
            m = 3;
            pm = [m,R1;
                m,R1,
                m,R1];
            I = bodyInertiaMatrix(pm);
            value = [[0,-z,y;
                    z,0,-x;
                    -y,x,0]*m*3,momentOfInertia(pm)];

            testCase.verifyLessThan(abs(I(4:6,:) - value), testCase.epsilon);
        end
    end
end

