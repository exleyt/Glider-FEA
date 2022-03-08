classdef testExcitingForces < matlab.unittest.TestCase

    properties
        epsilon
    end

    methods(TestMethodSetup)
        function defineEpsilon(testCase)
            testCase.epsilon = 1E-4;
        end
    end
    
    methods (Test)
        function testExcitingForce1(testCase)
            theta = pi/4;
            g = 9.81;
            k = 0.5;
            p = 1026; 
            FN = [0,0,0];
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,1,0;
                   0,0,1];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);
            value = [0,0,0,0,0,0];

            testCase.verifyLessThan(abs(F - value),testCase.epsilon)
        end

        function testExcitingForce2(testCase)
            theta = pi/4;
            g = 0;
            k = 0.5;
            p = 1026; 
            FN = [1,1,1]/norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,1,0;
                   0,0,1];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);
            value = [0,0,0,0,0,0];

            testCase.verifyLessThan(abs(F - value),testCase.epsilon)
        end

        function testExcitingForce3(testCase)
            theta = pi/4;
            g = 9.8;
            k = 0.5;
            p = 0; 
            FN = [1,1,1]/norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,1,0;
                   0,0,1];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);
            value = [0,0,0,0,0,0];

            testCase.verifyLessThan(abs(F - value),testCase.epsilon)
        end

        function testExcitingForce4(testCase)
            theta = pi/4;
            g = 9.8;
            k = 0;
            p = 1026; 
            FN = [2,5,4]/norm([2,5,4]);
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,0,0;
                   0,1,0];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);
            value = [FN(1)*0.5;
                     FN(2)*0.5;
                     FN(3)*0.5;
                     FN(3)/6;
                     -FN(3)/6;
                     (FN(2)-FN(1))/6]*p*g;

            testCase.verifyLessThan(abs(F - value),testCase.epsilon)
        end

        function testExcitingForce5(testCase)
            theta = pi/4;
            g = 9.8;
            k = 0;
            p = 1026; 
            FN = [1,1,1]/norm([1,1,1]);
            phi = [1;1;1;1;1;1];
            Tri = [1,0,0;
                   0,0,1;
                   0,1,0];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);
            value = [FN(1)*0.8660254037844389;
                     FN(2)*0.8660254037844389;
                     FN(3)*0.8660254037844389;
                     0;
                     0;
                     0]*p*g;

            testCase.verifyLessThan(abs(F(1:3) - value(1:3)),testCase.epsilon)
        end

        function testExcitingForce6(testCase)
            theta = pi/4;
            g = 9.8;
            k = 0;
            p = 1026; 
            FN = [2,4,3]/norm([2,4,3]);
            phi = [1;2;3;4;5;6];
            Tri = [1,0,0;
                   0,0,1;
                   0,1,0];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);

            r0 = Tri(1,:);
            ru = Tri(2,:) - r0;
            rv = Tri(3,:) - r0;

            s0 = k*(r0(3) - 1i*(r0(1)*cos(theta) + r0(2)*sin(theta)));
            su = k*(ru(3) - 1i*(ru(1)*cos(theta) + ru(2)*sin(theta)));
            sv = k*(rv(3) - 1i*(rv(1)*cos(theta) + rv(2)*sin(theta)));

            [L1,~,~] = surfIntPhiI(su,sv);
            value = [FN(1)*1i - phi(1)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(2)*1i - phi(2)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(3)*1i - phi(3)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     0;
                     0;
                     0]*p*g*-1i*(2*0.8660254037844389)*exp(s0)*L1;

            testCase.verifyLessThan(abs(F(1:3) - value(1:3)),testCase.epsilon)
        end

        function testExcitingForce7(testCase)
            theta = pi/4;
            g = 9.8;
            k = 0.6;
            p = 1026; 
            FN = [2,4,3]/norm([2,4,3]);
            phi = [1;2;3;4;5;6];
            Tri = [1,0,0;
                   0,0,1;
                   0,1,0];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);

            r0 = Tri(1,:);
            ru = Tri(2,:) - r0;
            rv = Tri(3,:) - r0;

            s0 = k*(r0(3) - 1i*(r0(1)*cos(theta) + r0(2)*sin(theta)));
            su = k*(ru(3) - 1i*(ru(1)*cos(theta) + ru(2)*sin(theta)));
            sv = k*(rv(3) - 1i*(rv(1)*cos(theta) + rv(2)*sin(theta)));

            [L1,~,~] = surfIntPhiI(su,sv);
            value = [FN(1)*1i - phi(1)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(2)*1i - phi(2)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(3)*1i - phi(3)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     0;
                     0;
                     0]*p*g*-1i*(2*0.8660254037844389)*exp(s0)*L1;

            testCase.verifyLessThan(abs(F(1:3) - value(1:3)),testCase.epsilon)
        end

        function testExcitingForce8(testCase)
            theta = pi/4;
            g = 9.8;
            k = 0.6;
            p = 1026; 
            FN = [2,4,3]/norm([2,4,3]);
            phi = [1;2;3;4;5;6];
            Tri = [1,0,0;
                   0,0,1;
                   0,1,0];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);

            r0 = Tri(1,:);
            ru = Tri(2,:) - r0;
            rv = Tri(3,:) - r0;

            s0 = k*(r0(3) - 1i*(r0(1)*cos(theta) + r0(2)*sin(theta)));
            su = k*(ru(3) - 1i*(ru(1)*cos(theta) + ru(2)*sin(theta)));
            sv = k*(rv(3) - 1i*(rv(1)*cos(theta) + rv(2)*sin(theta)));

            [L1,Lu,Lv] = surfIntPhiI(su,sv);
            value = [FN(1)*1i - phi(1)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(2)*1i - phi(2)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(3)*1i - phi(3)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     0;
                     0;
                     0]*p*g*-1i*(2*0.8660254037844389)*exp(s0)*L1;
            value = value + [0;
                     0;
                     0;
                     1i*(FN(3)*ru(2) - FN(2)*ru(3))*Lu + 1i*(FN(3)*rv(2) - FN(2)*rv(3))*Lv + ((FN(3)*r0(2) - FN(2)*r0(3)) - phi(4)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3)))*L1;
                     1i*(FN(1)*ru(3) - FN(3)*ru(1))*Lu + 1i*(FN(1)*rv(3) - FN(3)*rv(1))*Lv + ((FN(1)*r0(3) - FN(3)*r0(1)) - phi(5)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3)))*L1;
                     1i*(FN(2)*ru(1) - FN(1)*ru(2))*Lu + 1i*(FN(2)*rv(1) - FN(1)*rv(2))*Lv + ((FN(2)*r0(1) - FN(1)*r0(2)) - phi(6)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3)))*L1]*p*g*-1i*(2*0.8660254037844389)*exp(s0);
            testCase.verifyLessThan(abs(F(1:3) - value(1:3)),testCase.epsilon)
        end

        function testExcitingForce9(testCase)
            theta = pi/6;
            g = 9.8;
            k = 0.26;
            p = 1026; 
            FN = [2,6,3]/norm([2,6,3]);
            phi = [1;2;3;4;5;6];
            Tri = [1,0,0;
                   0,0,1;
                   0,1,0];
            [F] = excitingForce(Tri,phi,FN,p,k,g,theta);

            r0 = Tri(1,:);
            ru = Tri(2,:) - r0;
            rv = Tri(3,:) - r0;

            s0 = k*(r0(3) - 1i*(r0(1)*cos(theta) + r0(2)*sin(theta)));
            su = k*(ru(3) - 1i*(ru(1)*cos(theta) + ru(2)*sin(theta)));
            sv = k*(rv(3) - 1i*(rv(1)*cos(theta) + rv(2)*sin(theta)));

            [L1,Lu,Lv] = surfIntPhiI(su,sv);
            value = [FN(1)*1i - phi(1)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(2)*1i - phi(2)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     FN(3)*1i - phi(3)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3));
                     0;
                     0;
                     0]*p*g*-1i*(2*0.8660254037844389)*exp(s0)*L1;
            value = value + [0;
                     0;
                     0;
                     1i*(FN(3)*ru(2) - FN(2)*ru(3))*Lu + 1i*(FN(3)*rv(2) - FN(2)*rv(3))*Lv + ((FN(3)*r0(2) - FN(2)*r0(3)) - phi(4)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3)))*L1;
                     1i*(FN(1)*ru(3) - FN(3)*ru(1))*Lu + 1i*(FN(1)*rv(3) - FN(3)*rv(1))*Lv + ((FN(1)*r0(3) - FN(3)*r0(1)) - phi(5)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3)))*L1;
                     1i*(FN(2)*ru(1) - FN(1)*ru(2))*Lu + 1i*(FN(2)*rv(1) - FN(1)*rv(2))*Lv + ((FN(2)*r0(1) - FN(1)*r0(2)) - phi(6)*k*(FN(1)*cos(theta) + FN(2)*sin(theta) + 1i*FN(3)))*L1]*p*g*-1i*(2*0.8660254037844389)*exp(s0);
            testCase.verifyLessThan(abs(F(1:3) - value(1:3)),testCase.epsilon)
        end
    end
end

