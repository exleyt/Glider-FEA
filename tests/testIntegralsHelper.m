classdef testIntegralsHelper  
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
    end

    methods (Access = public)

        function [b1,bu,bv] = evaluateIntegral(self,su,sv)
            bv1 = @(v) exp(sv*v);
            bu1 = @(u) exp(su*u)*self.pasq(bv1,0,1 - u);
            b1 = self.pasq(bu1,0,1);

            buu = @(u) u*exp(su*u)*self.pasq(bv1,0,1 - u);
            bu = self.pasq(buu,0,1);

            bvv = @(v) v*exp(sv*v);
            bvu = @(u) exp(su*u)*self.pasq(bvv,0,1 - u);
            bv = self.pasq(bvu,0,1);
        end
    end
end

