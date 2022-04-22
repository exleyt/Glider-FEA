classdef greenComparisonTests
    
    methods (Access = private)
        function [result] = pasq(self,R,sx3,K,a,b,e,dmax)
            fa = besselj(0,a*R)*exp(a*sx3)/(a-K);
            fb = besselj(0,b*R)*exp(b*sx3)/(b-K);
            [m,fm,whole] = self.pasqm(R,sx3,K,a,fa,b,fb);
            result = self.pasqr(R,sx3,K,a,fa,m,fm,b,fb,whole,e,0,dmax);
        end 

        function [m,fm,whole] = pasqm(~,R,sx3,K,a,fa,b,fb)
            m = (a + b)/2;
            fm = besselj(0,m*R)*exp(m*sx3)/(m-K);
            whole = abs(b - a)*(fa + 4*fm + fb)/6;
        end

        function [result] = pasqr(self,R,sx3,K,a,fa,m,fm,b,fb,whole,e,d,dmax)
            [lm,flm,left] = self.pasqm(R,sx3,K,a,fa,m,fm);
            [rm,frm,right] = self.pasqm(R,sx3,K,m,fm,b,fb);
            delta = left + right - whole;
            if abs(delta) <= 15 * e
                result = left + right + delta / 15;
            elseif d >= dmax
                result = left + right + delta / 15;
            else
                result = self.pasqr(R,sx3,K,a,fa,lm,flm,m,fm,left,e/2,d+1,dmax);
                result = result + self.pasqr(R,sx3,K,m,fm,rm,frm,b,fb,right,e/2,d+1,dmax); 
            end
        end

        function [result] = green(self,x,xi,K,es,ep,infi,dmax)
            dx2 = (x - xi).^2;
            r1 = 1/sqrt(dx2(1) + dx2(2) + dx2(3));
            sx3 = x(3) + xi(3);
            r2 = 1/sqrt(dx2(1) + dx2(2) + sx3^2);
            R = sqrt(dx2(1) + dx2(2));
            result = r1 + r2 + 2*K*(self.pasq(R,sx3,K,0,K-es,ep,dmax) + ...
                self.pasq(R,sx3,K,K+es,infi,ep,dmax)) - ...
                2*pi*1i*K*exp(K*sx3)*besselj(0,K*R);
        end

        function [result] = greenPartialDot(self,x,xi,K,j,f,e,es,ep,infi,dmax)
            xi(j) = xi(j) + e;
            fe = self.green(x,xi,K,es,ep,infi,dmax);
            result = (fe - f) / e;
        end

        function [f,df] = greenAndPartialDot(self,x,xi,FN,K,e,es,ep,infi,dmax)
            f = self.green(x,xi,K,es,ep,infi,dmax); % the green function evaluated at xi
            df = 0;
            for j = 1:3
                df = df + FN(j)*self.greenPartialDot(x,xi,K,j,f,e,es,ep,infi,dmax);
            end
        end
        
        function [f,df] = greenAndPartial(self,x,xi,FN,K,e,es,ep,infi,dmax)
            f = self.green(x,xi,K,es,ep,infi,dmax); % the green function evaluated at xi 
            fe = self.green(x,xi + e*FN,K,es,ep,infi,dmax);
            df = (fe - f) / e; 
        end

        function [x,xi,FN,K] = randGreenInputs(~)
            x = [rand(1),rand(1),-abs(rand(1))];
            xi = [rand(1),rand(1),-abs(rand(1))];
            FN = [rand(1),rand(1),rand(1)];
            FN = FN / norm(FN);
            K = rand(1);
        end
    end

    methods (Access = public)
        % They converge around E = [1E-10,1E-4]
        function runPartialVSDot(self)
            [x,xi,FN,K] = self.randGreenInputs();

            dmax = 20;
            infi = 100;
            es = 1E-10;

            EP = [1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1];
            E = [1E-14,1E-13,1E-12,1E-11,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1];

            
            for ep = EP
                for e = E
                    [~,df1] = self.greenAndPartialDot(x,xi,FN,K,e,es,ep,infi,dmax);
                    [~,df2] = self.greenAndPartial(x,xi,FN,K,e,es,ep,infi,dmax);

                    avg = (df1 + df2)/2;
                    dev = abs(real(df1 - avg)) + 1i*abs(imag(df1 - avg));
                    pdev = real(dev)/abs(real(avg)) + 1i*imag(dev)/abs(imag(avg));

                    fprintf(1,"%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",ep,e,real(pdev),imag(pdev),real(avg),imag(avg));
                end
            end
        end

        function runPartial(self)
            [x,xi,FN,K] = self.randGreenInputs();
            e = 1E-6; % A consistent value from runPartialVSDot
            EP = [1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1];
            ES = [1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1];
            INFI = [5,1E2,50,1E3,500,1E4];
            dmax = 20;

            for ep = EP
                for es = ES
                    for infi = INFI
                        [f,df] = self.greenAndPartial(x,xi,FN,K,e,es,ep,infi,dmax);

                        fprintf(1,"%.11f,%.11f,%d,%.11f,%.11f,%.11f\n",ep,es,infi,real(f),real(df),imag(df));
                    end
                end
            end
        end

        function runPartialDepth(self,e,ep,es,infi)
            [x,xi,FN,K] = self.randGreenInputs();
            for dmax = [5,10,15,20,25,30,35,40,45,50]
                [f,df] = self.greenAndPartial(x,xi,FN,K,e,es,ep,infi,dmax);
                fprintf(1,"%d,%.11f,%.11f,%.11f,%.11f\n",dmax,real(f),imag(f),real(df),imag(df));
            end
        end
    end
end

