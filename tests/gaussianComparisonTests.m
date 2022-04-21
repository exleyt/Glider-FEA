classdef gaussianComparisonTests    

    methods (Access = public)
        function [result] = pasq(self,f,a,b)
            fa = f(a);
            fb = f(b);
            [m,fm,whole] = self.pasqm(f,a,fa,b,fb);
            e = 1E-2;
            result = self.pasqr(f,a,fa,m,fm,b,fb,whole,e,0);
        end 

        function [result] = pasq2(self,f,u,a,b)
            fa = f(u,a);
            fb = f(u,b);
            [m,fm,whole] = self.pasqm2(f,u,a,fa,b,fb);
            e = 1E-2;
            result = self.pasqr2(f,u,a,fa,m,fm,b,fb,whole,e,0);
        end

        function [m,fm,whole] = pasqm(~,f,a,fa,b,fb)
            m = (a + b)/2;
            fm = f(m);
            whole = abs(b - a)*(fa + 4*fm + fb)/6;
        end

        function [m,fm,whole] = pasqm2(~,f,u,a,fa,b,fb)
            m = (a + b)/2;
            fm = f(u,m);
            whole = abs(b - a)*(fa + 4*fm + fb)/6;
        end

        function [result] = pasqr(self,f,a,fa,m,fm,b,fb,whole,e,d)
            [lm,flm,left] = self.pasqm(f,a,fa,m,fm);
            [rm,frm,right] = self.pasqm(f,m,fm,b,fb);
            delta = left + right - whole;
            if abs(delta) <= 15 * e
                result = left + right + delta / 15;
            elseif d >= 20
                result = left + right + delta / 15;
            else
                result = self.pasqr(f,a,fa,lm,flm,m,fm,left,e/2,d+1);
                result = result + self.pasqr(f,m,fm,rm,frm,b,fb,right,e/2,d+1); 
            end
        end

        function [result] = pasqr2(self,f,u,a,fa,m,fm,b,fb,whole,e,d)
            [lm,flm,left] = self.pasqm2(f,u,a,fa,m,fm);
            [rm,frm,right] = self.pasqm2(f,u,m,fm,b,fb);
            delta = left + right - whole;
            if abs(delta) <= 15 * e
                result = left + right + delta / 15;
            elseif d >= 20
                result = left + right + delta / 15;
            else
                result = self.pasqr2(f,u,a,fa,lm,flm,m,fm,left,e/2,d+1);
                result = result + self.pasqr2(f,u,m,fm,rm,frm,b,fb,right,e/2,d+1); 
            end
        end
    
        function [df] = getGreenPartialXINormal(~,x,xi,FN,K)
            [~,df] = greenFunctionAndPartialXINormal(x,xi,FN,K);
        end

        function [result] = getGreen0(~,x,xi)
            dx2 = (x - xi).^2;
            r1 = 1/sqrt(dx2(1) + dx2(2) + dx2(3));
            sx3 = x(3) + xi(3);
            r2 = 1/sqrt(dx2(1) + dx2(2) + sx3^2);
            result = r1 + r2;
        end

        function [f,df] = getGreenAndPartialXINormal0(self,x,xi,FN)
            f = self.getGreen0(x,xi);
            e = 10^-20; % epsilon for estimating 
            fe = self.getGreen0(x,xi + e*FN);
            df = (fe - f)/e; 
        end

        function [df] = getGreenPartialXINormal0(self,x,xi,FN)
            [~,df] = self.getGreenAndPartialXINormal0(x,xi,FN);
        end
    end

    methods (Access = private)
        function [T,CP,FN] = getMesh(self,path)
            mesh = stlread(path);
            to = triangulation(mesh.ConnectivityList, ...
                mesh.Points*1E4); %+ [.5,0,0]);
            
            [T,CP,FN] = self.workMesh(to);
        end

        function [T,CP,FN] = workMesh(~,mesh)            
            [S,~] = size(mesh.ConnectivityList);
            
            T = zeros(3,3,S);
            for j = 1:S
                T(:,:,j) = mesh.Points(mesh.ConnectivityList(j,:),:);
            end
            
            CP = incenter(mesh);
            FN = -faceNormal(mesh);
        end

        function [G] = getGreenEstimateG(~,T,CP,f,w,s,K)
            [~,~,N] = size(T);
            G = zeros(N,N);

            for k = 1:N
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = f(:,1)*ru + f(:,2)*rv + r0;
                A = 0.5*norm(cross(ru,rv));

                for n = 1:N
                    if n ~= k
                        for m = 1:s
                            [fg] = greenFunction(CP(n,:),r(m,:),K);
                            G(n,k) = G(n,k) + fg*w(m);
                        end
                        G(n,k) = G(n,k)*A;
                    end
                end
            end
        end

        function [G] = getGreenRealG(self,T,CP,K)
            [~,~,N] = size(T);
            G = zeros(N,N);
                                     
            parfor k = 1:N
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = @(u,v) u*ru + v*rv + r0;
                A = norm(cross(ru,rv));

                for n = 1:N
                    if n ~= k
                        Gv = @(u,v) greenFunction(CP(n,:),r(u,v),K);
                        Gu = @(u) self.pasq2(Gv,u,0,1-u);
                        G(n,k) = A*self.pasq(Gu,0,1);
                    end
                end
            end
        end

        function [M] = getGreenRealM(self,T,CP,FN,K)
            [~,~,N] = size(T);
            M = zeros(N,N);
                                     
            parfor k = 1:N
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = @(u,v) u*ru + v*rv + r0;
                A = norm(cross(ru,rv));
                for n = 1:N
                    if n ~= k
                        Gv = @(u,v) self.getGreenPartialXINormal(CP(n,:),r(u,v),FN(k,:),K);
                        Gu = @(u) self.pasq2(Gv,u,0,1-u);
                        M(n,k) = A*self.pasq(Gu,0,1);
                    else
                        M(n,k) = 2*pi;
                    end
                end
            end
        end

        function [G,M] = getGreenEstimate(~,T,CP,FN,f,w,s,K)
            [~,~,N] = size(T);
            G = zeros(N,N);
            M = zeros(N,N);

            for k = 1:N
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = f(:,1)*ru + f(:,2)*rv + r0;
                A = 0.5*norm(cross(ru,rv));

                for n = 1:N
                    if k ~= n
                        for m = 1:s
                            [fg,fm] = greenFunctionAndPartialXINormal(CP(n,:),r(m,:),FN(k,:),K);
                            G(n,k) = G(n,k) + fg*w(m);
                            M(n,k) = M(n,k) + fm*w(m);
                        end
                        G(n,k) = G(n,k)*A;
                        M(n,k) = M(n,k)*A;
                    else
                        M(n,k) = 2*pi;
                    end
                end
            end
        end

        function [f,w,s] = getGaussians(~)
            f = zeros(20,2,12);
            w = zeros(20,12);
            
            f(:,:,1) = [1/6,1/6;
                        2/3,1/6;
                        1/6,2/3;
                        zeros(17,2)];
            w(:,1) =   [1/3;
                        1/3;
                        1/3;
                        zeros(17,1)];
            s(1) = 3;
            
            f(:,:,2) = [0.1012865073235,0.1012865073235; 
                        0.7974269853531,0.1012865073235;
                        0.1012865073235,0.7974269853531;
                        0.4701420641051,0.0597158717898;
                        0.4701420641051,0.4701420641051;
                        0.0597158717898,0.4701420641051;
                        0.3333333333333,0.3333333333333;
                        zeros(13,2)];
            w(:,2) =   [0.1259391805448;
                        0.1259391805448;
                        0.1259391805448;
                        0.1323941527885;
                        0.1323941527885;
                        0.1323941527885;
                        0.225;
                        zeros(13,1)];
            s(2) = 7;
            
            f(:,:,3) = [0.0651301029022,0.0651301029022;
                        0.8697397941956,0.0651301029022;
                        0.0651301029022,0.8697397941956;
                        0.3128654960049,0.0486903154253;
                        0.6384441885698,0.3128654960049;
                        0.0486903154253,0.6384441885698;
                        0.6384441885698,0.0486903154253;
                        0.3128654960049,0.6384441885698;
                        0.0486903154253,0.0486903154253;
                        0.2603459660790,0.2603459660790;
                        0.4793080678419,0.2603459660790;
                        0.2603459660790,0.4793080678419;
                        1/3,1/3;
                        zeros(7,2)];
            w(:,3) =   [0.0533472356088;
                        0.0533472356088;
                        0.0533472356088;
                        0.0771137608903;
                        0.0771137608903;
                        0.0771137608903;
                        0.0771137608903;
                        0.0771137608903;
                        0.0771137608903;
                        0.1756152574332;
                        0.1756152574332;
                        0.1756152574332;
                        -0.1495700444677;
                        zeros(7,1)];
            s(3) = 13;
            
            f(:,:,4) = [1/3,1/3;
                        zeros(19,2)];
            w(:,4) =   [1;
                        zeros(19,1)];
            s(4) = 1;
            
            f(:,:,5) = [0,0.5;
                        0.5,0;
                        0.5,0.5;
                        zeros(17,2)];
            w(:,5) =   [1/3;
                        1/3;
                        1/3;
                        zeros(17,1)];
            s(5) = 3;
            
            f(:,:,6) = [1/3,1/3;
                        1/5,1/5;
                        1/5,3/5;
                        3/5,1/5;
                        zeros(16,2)];
            w(:,6) =  [-27/48;
                        25/48;
                        25/48;
                        25/48;
                        zeros(16,1)];
            s(6) = 4;
            
            f(:,:,7) = [1/3,1/3;
                        2/15,11/15;
                        2/15,2/15;
                        11/15,2/15;
                        zeros(16,2)];
            w(:,7) =  [-27/48;
                        25/48;
                        25/48;
                        25/48;
                        zeros(16,1)];
            s(7) = 4;
            
            f(:,:,8) = [0.445948490915970,0.445948490915970;
                        0.445948490915970,0.108103018168070;
                        0.108103018168070,0.445948490915970;
                        0.091576213509770,0.091576213509770;
                        0.091576213509770,0.816847572980460;
                        0.816847572980460,0.091576213509770;
                        zeros(14,2)];
            w(:,8) =   [0.223381589678010;
                        0.223381589678010;
                        0.223381589678010;
                        0.109951743655320;
                        0.109951743655320;
                        0.109951743655320;
                        zeros(14,1)];
            s(8) = 6;
            
            f(:,:,9) = [0.333333333333330,0.333333333333330;
                        0.470142064105110,0.470142064105110;
                        0.470142064105110,0.059715871789770;
                        0.059715871789770,0.470142064105110;
                        0.101286507323460,0.101286507323460;
                        0.101286507323460,0.797426985353090;
                        0.797426985353090,0.101286507323460;
                        zeros(13,2)];
            w(:,9) =   [0.225000000000000;
                        0.132394152788510;
                        0.132394152788510;
                        0.132394152788510;
                        0.125939180544830;
                        0.125939180544830;
                        0.125939180544830;
                        zeros(13,1)];
            s(9) = 7;
            
            f(:,:,10) = [0.249286745170910,0.249286745170910;
                        0.249286745170910,0.501426509658180;
                        0.501426509658180,0.249286745170910;
                        0.063089014491500,0.063089014491500;
                        0.063089014491500,0.873821971017000;
                        0.873821971017000,0.063089014491500;
                        0.310352451033780,0.636502499121400;
                        0.636502499121400,0.053145049844820;
                        0.053145049844820,0.310352451033780;
                        0.636502499121400,0.310352451033780;
                        0.310352451033780,0.053145049844820;
                        0.053145049844820,0.636502499121400;
                        zeros(8,2)];
            w(:,10) =   [0.116786275726380;
                        0.116786275726380;
                        0.116786275726380;
                        0.050844906370210;
                        0.050844906370210;
                        0.050844906370210;
                        0.082851075618370;
                        0.082851075618370;
                        0.082851075618370;
                        0.082851075618370;
                        0.082851075618370;
                        0.082851075618370;
                        zeros(8,1)];
            s(10) = 12;
            
            f(:,:,11) = [0.33333333333333,0.333333333333330;
                        0.260345966079040,0.260345966079040;
                        0.260345966079040,0.479308067841920;
                        0.479308067841920,0.260345966079040;
                        0.065130102902220,0.065130102902220;
                        0.065130102902220,0.869739794195570;
                        0.869739794195570,0.065130102902220;
                        0.312865496004870,0.638444188569810;
                        0.638444188569810,0.048690315425320;
                        0.048690315425320,0.312865496004870;
                        0.638444188569810,0.312865496004870;
                        0.312865496004870,0.048690315425320;
                        0.048690315425320,0.638444188569810;
                        zeros(7,2)];
            w(:,11) =  [-0.14957004446768;
                        0.175615257433210;
                        0.175615257433210;
                        0.175615257433210;
                        0.053347235608840;
                        0.053347235608840;
                        0.053347235608840;
                        0.077113760890260;
                        0.077113760890260;
                        0.077113760890260;
                        0.077113760890260;
                        0.077113760890260;
                        0.077113760890260;
                        zeros(7,1)];
            s(11) = 13;
            
            f(:,:,12) = [0.33333333333333,0.333333333333330;
                        0.459292588292720,0.459292588292720;
                        0.459292588292720,0.081414823414550;
                        0.081414823414550,0.459292588292720;
                        0.170569307751760,0.170569307751760;
                        0.170569307751760,0.658861384496480;
                        0.658861384496480,0.170569307751760;
                        0.050547228317030,0.050547228317030;
                        0.050547228317030,0.898905543365940;
                        0.898905543365940,0.050547228317030;
                        0.263112829634640,0.728492392955400;
                        0.728492392955400,0.008394777409960;
                        0.008394777409960,0.263112829634640;
                        0.728492392955400,0.263112829634640;
                        0.263112829634640,0.008394777409960;
                        0.008394777409960,0.728492392955400;
                        zeros(4,2)];
            w(:,12) =  [0.144315607677790;
                        0.095091634267280;
                        0.095091634267280;
                        0.095091634267280;
                        0.103217370534720;
                        0.103217370534720;
                        0.103217370534720;
                        0.032458497623200;
                        0.032458497623200;
                        0.032458497623200;
                        0.027230314174430;
                        0.027230314174430;
                        0.027230314174430;
                        0.027230314174430;
                        0.027230314174430;
                        0.027230314174430;
                        zeros(4,1)];
            s(12) = 16;
        end

        function [phi] = getVelocityPotential(~,G,M,FN6)
            [N,~] = size(G);
            Gsum = zeros(N,6);

            for m = 1:N
                for n = 1:N
                    Gsum(m,:) = Gsum(m,:) + G(m,n)*FN6(m,:);                     
                end
            end

            phi = zeros(N,6);
            for j = 1:6
                phi(:,j) = linsolve(M,Gsum(:,j));
            end
        end

        function [G,M] = getGreenReal0(self,T,CP,FN)
            [~,~,N] = size(T);
            M = zeros(N,N);
            G = zeros(N,N);
                                     
            parfor k = 1:N
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = @(u,v) u*ru + v*rv + r0;
                A = norm(cross(ru,rv));
                for n = 1:N
                    if n ~= k
                        Mv = @(u,v) self.getGreenPartialXINormal0(CP(n,:),r(u,v),FN(k,:));
                        Mu = @(u) self.pasq2(Mv,u,0,1-u);
                        M(n,k) = A*self.pasq(Mu,0,1);

                        Gv = @(u,v) self.getGreen0(CP(n,:),r(u,v));
                        Gu = @(u) self.pasq2(Gv,u,0,1-u);
                        G(n,k) = A*self.pasq(Gu,0,1);
                    else
                        M(n,k) = 2*pi;
                    end
                end
            end
        end

        function [G,M] = getGreenEstimate0(self,T,CP,FN,f,w,s)
            [~,~,N] = size(T);
            G = zeros(N,N);
            M = zeros(N,N);

            parfor k = 1:N
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = f(:,1)*ru + f(:,2)*rv + r0;
                A = 0.5*norm(cross(ru,rv));

                for n = 1:N
                    if k ~= n
                        for m = 1:s
                            [fg,fm] =   self.getGreenAndPartialXINormal0(CP(n,:),r(m,:),FN(k,:));
                            G(n,k) = G(n,k) + fg*w(m);
                            M(n,k) = M(n,k) + fm*w(m);
                        end
                        G(n,k) = G(n,k)*A;
                        M(n,k) = M(n,k)*A;
                    else
                        M(n,k) = 2*pi;
                    end
                end
            end
        end
    end

    methods (Access = public)
        % Runs the tests for the selected Gaussians where select is a list
        %  of the indecies of s to be run
        function runSelectedG(self,select)
            model = createpde();
            importGeometry(model,"test_models\sphere geometry.stl");
            scale(model.Geometry, 0.001);
            generateMesh(model, 'GeometricOrder','linear','Hmin',0.01,'Hmax',0.1);
            P = model.Mesh.Nodes(:,model.Mesh.Elements(:,1)).';
            CL = [1,2,3;2,3,4;3,4,1;4,1,2];
            
            [T,CP,~] = self.workMesh(triangulation(CL,P));
            [~,~,N] = size(T);
            [f,w,s] = self.getGaussians();

            K = [1.0071,0.2518,0.0403];

            for k = K
                GR = self.getGreenRealG(T,CP,k);

                for j = select
                    G = self.getGreenEstimateG(T,CP,f(:,:,j),w(:,j),s(j),k);
                    
                    diff = (G - GR)./GR;
                    for i = 1:N
                        diff(i,i) = 0;
                    end
                    diff = abs(real(diff)) + 1i*abs(imag(diff));

                    avg = sum(sum(diff))/(N*N - N); % Not counting the diags
                    maxi = max(max(diff));
                    mini = min(min(diff + diag(zeros(1,N) + 100)));

                    fprintf(1,"%d,%f,%f%+fi,%f%+fi,%f%+fi\n",j,k,real(avg),imag(avg),real(mini),imag(mini),real(maxi),imag(maxi));
                end
            end
        end

        function runAllG(self)
            [~,~,s] = self.getGaussians();
            [~,S] = size(s);
            self.runSelectedG(1:S);
        end

        function runBestsG(self)
            self.runSelectedG([1,5,8,10]);
        end

        function runSelectedM(self,select)
            model = createpde();
            importGeometry(model,"test_models\sphere geometry.stl");
            scale(model.Geometry, 0.001);
            generateMesh(model, 'GeometricOrder','linear','Hmin',0.01,'Hmax',0.1);
            P = model.Mesh.Nodes(:,model.Mesh.Elements(:,1)).';
            CL = [1,2,3;2,3,4;3,4,1;4,1,2];
            
            [T,CP,FN] = self.workMesh(triangulation(CL,P));
            [~,~,N] = size(T);
            [f,w,s] = self.getGaussians();

            K = [1.0071,0.2518,0.0403];

            for k = K
                MR = self.getGreenRealM(T,CP,FN,k);

                for j = select
                    [~,M] = self.getGreenEstimate(T,CP,FN,f(:,:,j),w(:,j),s(j),k);
                    
                    diff = (M - MR)./MR;
                    for i = 1:N
                        diff(i,i) = 0;
                    end
                    diff = abs(real(diff)) + 1i*abs(imag(diff));

                    avg = sum(sum(diff))/(N*N - N); % Not counting the diags
                    maxi = max(max(diff));
                    mini = min(min(diff + diag(zeros(1,N) + 100)));

                    fprintf(1,"%d,%f,%f%+fi,%f%+fi,%f%+fi\n",j,k,real(avg),imag(avg),real(mini),imag(mini),real(maxi),imag(maxi));
                end
            end
        end

        function runAllM(self)
            [~,~,s] = self.getGaussians();
            [~,S] = size(s);
            self.runSelectedM(1:S);
        end

        function runBestsM(self)
            self.runSelectedM([1,5,8,10]);
        end

        function runSelected0(self,select)
            model = createpde();
            importGeometry(model,"test_models\sphere geometry.stl");
            scale(model.Geometry, 0.001);
            generateMesh(model, 'GeometricOrder','linear','Hmin',0.01,'Hmax',0.1);
            to = triangulation(model.Mesh.Elements.',model.Mesh.Nodes.');
            [CL,P] = freeBoundary(to);

            [T,CP,FN] = self.workMesh(triangulation(CL,P));
            FN6 = normal6DOF(CP,FN);
            [~,~,N] = size(T);
            [f,w,s] = self.getGaussians();

            [GR,MR] = self.getGreenReal0(T,CP,FN);            
            phiR = self.getVelocityPotential(GR,MR,FN6);

            for j = select
                [G,M] = self.getGreenEstimate0(T,CP,FN,f(:,:,j),w(:,j),s(j));
                phi = self.getVelocityPotential(G,M,FN6); % (N,6)

                diff = (phiR - phi)./(phiR + 1E-20);
                diff = abs(real(diff)) + 1i*abs(imag(diff));
                avg = sum(sum(diff))/(N*6);
                maxi = max(max(diff));
                mini = min(min(diff));
                
                fprintf(1,"%d,%f%+fi,%f%+fi,%f%+fi\n",j,real(avg),imag(avg),real(mini),imag(mini),real(maxi),imag(maxi));
            end
        end

        function runAll0(self)
            [~,~,s] = self.getGaussians();
            [~,S] = size(s);
            self.runSelected0(1:S);
        end

        function runBests0(self)
            self.runSelected0([1,5,8,10]);
        end
    end
end