classdef guassianComparisonTests    

    methods (Access = private)
         function [T,CP] = getMesh(self,path)
            mesh = stlread(path);
            to = triangulation(mesh.ConnectivityList, ...
                mesh.Points*1E4); %+ [.5,0,0]);

            [S,~] = size(mesh.ConnectivityList);
            
            T = zeros(3,3,S);
            for j = 1:S
                T(:,:,j) = to.Points(to.ConnectivityList(j,:),:);
            end

            CP  = incenter(to);
         end
        
         function [G] = getGreenEstimate0Diag(self,T,CP,f,w,s,K)
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

         function [G] = getGreenEstimateDiag(self,T,CP,f,w,s,K)
            [~,~,N] = size(T);
            G = zeros(N,1);

            for k = 1:N
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = f(:,1)*ru + f(:,2)*rv + r0;
                A = 0.5*norm(cross(ru,rv));

                for m = 1:s
                    [fg] = greenFunction(CP(k,:),r(m,:),K);
                    G(k) = G(k) + fg*w(m);
                end
                G(k) = G(k)*A;
            end
        end

        function [GR] = getGreenReal0Diag(self,T,CP,K,helper)
            [~,~,N] = size(T);
            GR = zeros(N,N);
                                     
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
                        Gu = @(u) helper.pasq2(Gv,u,0,1-u);
                        GR = A*helper.pasq(Gu,0,1);
                    end
                end
            end
        end

        function [GR] = getGreenRealDiag(self,T,CP,K,helper)
            [~,~,N] = size(T);
            GR = zeros(N,1);
                                     
            for k = 1:N
                disp(k)
                Tk = T(:,:,k);
                r0 = Tk(1,:);
                ru = Tk(2,:) - r0;
                rv = Tk(3,:) - r0;
                r = @(u,v) u*ru + v*rv + r0;
                A = norm(cross(ru,rv));

                 Gv = @(u,v) greenFunction(CP(k,:),r(u,v),K);
                Gu = @(u) helper.pasq2(Gv,u,0,1-u);
                GR = A*helper.pasq(Gu,0,1);
            end
        end

        function [f,w,s] = getGuassians(self)
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
    end

    methods (Access = public)
        function run(self)
            [T,CP] = self.getMesh("test_models\sphere mesh.stl");
            [~,~,N] = size(T);
            [f,w,s] = self.getGuassians();
            [~,S] = size(s);
            helper = testIntegralsHelper;

            for K = 0.4
                GRD = self.getGreenRealDiag(T,CP,K,helper);
                GR0D = self.getGreenReal0Diag(T,CP,K,helper);

                for j = 1:S
                    GD = self.getGreenEstimateDiag(T,CP,f(:,:,j),w(:,j),s(j),K);
                    G0D = self.getGreenEstimate0Diag(T,CP,f(:,:,j),w(:,j),s(j),K);
                    
                    diff = abs(GD - GRD);
                    diff0 = abs(G0D - GR0D);

                    avg = sum(sum(diff))/(N*N);
                    maxi = max(max(diff));
                    mini = min(min(diff));

                    avg0 = sum(sum(diff0))/(N*N);
                    maxi0 = max(max(diff0));
                    mini0 = min(min(diff0));

                    disp([j,avg,mini,maxi])
                    disp([j,avg0,mini0,maxi0])
                end
            end
        end
    end
end

