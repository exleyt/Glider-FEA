function [v,s,t,p,rt] = asq_test(fun,a,v,e,dmax)
% Compares the output of asq to the builtin matlab integral
%
% It does this by integrating fun from a to b with a maximum asq error 
% of e and maximum asq depth dmax
%
% It then returns if the asq value match the int values in a boolean v, 
% informative testing text in string s, if the performance time of asq is 
% better than int in a boolean t, a performance string that tells the 
% time difference between the two calculations in p, and the time value
% iteself is in rt
    r1 = @() adaptive_simpson_quadrature(fun,a,v,e,dmax);
    r2 = @() integral(fun,a,v);
    rt = timeit(r2) - timeit(r1);
    t = rt >= 0;
    p = ['ASQ Performance:', char(fun),'[',sprintf('%.2f',a),',',sprintf('%.2f',v),']',': ', sprintf('%.4f',rt)];
    if (~t)
        disp(p);
    end
    if r1() <= r2() + e && r1() >= r2() - e
        s = ['ASQ Test:', char(fun),'[',sprintf('%.2f',a),',',sprintf('%.2f',v),']',': Passed'];
        v = 1;
    else
        s = ['ASQ Test:', char(fun),'[',sprintf('%.2f',a),',',sprintf('%.2f',v),']',': Failed'];
        disp(s);
        v = 0;
    end
end

